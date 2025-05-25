#!/usr/bin/env python3

import argparse
import boto3
import os
import re
import sys
import threading
import time
from botocore.exceptions import ClientError, ReadTimeoutError
from concurrent.futures import ThreadPoolExecutor, as_completed
import pickle
import hashlib

# FAISS and embedding imports
try:
    import faiss
except ImportError:
    print(f"Error: FAISS not installed. Run: pip install faiss-gpu-cu12 sentence-transformers pymupdf", file=sys.stderr)
    print(f"For CPU-only FAISS: pip install faiss-cpu sentence-transformers pymupdf", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error: FAISS installation appears corrupted. Try reinstalling:", file=sys.stderr)
    print(f"  pip uninstall faiss-gpu-cu12 faiss-cpu", file=sys.stderr)
    print(f"  pip install faiss-cpu sentence-transformers pymupdf", file=sys.stderr)
    print(f"Original error: {e}", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
    from sentence_transformers import SentenceTransformer
    import fitz  # PyMuPDF
except ImportError as e:
    print(f"Error: Required packages not installed. Run: pip install pymupdf sentence-transformers", file=sys.stderr)
    sys.exit(1)

# --- Configuration Constants ---
# Default Model ID - will be auto-selected based on context size
DEFAULT_BEDROCK_MODEL_ID = "us.amazon.nova-micro-v1:0"
# Default AWS Region for Bedrock and S3 services (will be overridden by profile region if available)
DEFAULT_AWS_REGION = "us-east-1"

# Document processing limits - now used for chunking strategy
MAX_CHUNK_SIZE = 3000  # Characters per chunk for vector search (increased for fewer chunks)
CHUNK_OVERLAP = 500    # Overlap between chunks (proportionally increased)
TOP_K_CHUNKS = None    # Number of most relevant chunks to send to LLM (None = all chunks)

# Embedding cache settings
ENABLE_EMBEDDING_CACHE = True
CACHE_DIR = ".embedding_cache"

# Bedrock inference parameters and model limits
MODEL_LIMITS = {
    "us.amazon.nova-micro-v1:0": 128000,    # 128k tokens
    "us.amazon.nova-lite-v1:0": 300000,     # 300k tokens  
    "us.amazon.nova-premier-v1:0": 1000000, # 1M tokens
    "us.amazon.nova-pro-v1:0": 300000       # 300k tokens (fallback)
}

DEFAULT_INPUT_TOKENS = 30000  # Placeholder for typical model context window
DEFAULT_OUTPUT_TOKENS = 10000  # Default output token limit for model responses
DEFAULT_TOP_P = 0.9
DEFAULT_TEMPERATURE = 0.2


class Spinner:
    """Simple spinner for showing progress during long operations"""
    def __init__(self, message="Processing"):
        self.message = message
        self.spinning = False
        self.spinner_thread = None
        
    def start(self):
        """Start the spinner"""
        self.spinning = True
        self.spinner_thread = threading.Thread(target=self._spin)
        self.spinner_thread.daemon = True
        self.spinner_thread.start()
        
    def stop(self):
        """Stop the spinner"""
        self.spinning = False
        if self.spinner_thread:
            self.spinner_thread.join()
        # Clear the spinner line
        print("\r" + " " * (len(self.message) + 10) + "\r", end="", flush=True)
        
    def _spin(self):
        """Internal method to handle the spinning animation"""
        chars = "|/-\\"
        idx = 0
        while self.spinning:
            print(f"\r{self.message} {chars[idx % len(chars)]}", end="", flush=True)
            idx += 1
            time.sleep(0.1)


def detect_gpu():
    """Detect if NVIDIA GPU is available for FAISS and if FAISS has GPU support"""
    # First check if FAISS has GPU support
    try:
        faiss.StandardGpuResources()
        # If we get here, FAISS has GPU support
        gpu_available = False
        
        # Check if GPU hardware is available
        try:
            import torch
            if torch.cuda.is_available():
                gpu_count = torch.cuda.device_count()
                gpu_name = torch.cuda.get_device_name(0)
                print(f"GPU detected: {gpu_name} (using GPU-accelerated FAISS)")
                gpu_available = True
        except ImportError:
            pass
        
        # Alternative check using nvidia-ml-py
        if not gpu_available:
            try:
                import pynvml
                pynvml.nvmlInit()
                gpu_count = pynvml.nvmlDeviceGetCount()
                if gpu_count > 0:
                    print(f"GPU detected: {gpu_count} NVIDIA GPU(s) (using GPU-accelerated FAISS)")
                    gpu_available = True
            except ImportError:
                pass
        
        if not gpu_available:
            print("GPU hardware not available, using CPU FAISS")
        
        return gpu_available
        
    except (AttributeError, ImportError):
        # FAISS doesn't have GPU support (faiss-cpu installed)
        # Don't print here - will print only if FAISS is actually used
        return False


def extract_text_from_markdown(md_path):
    """Extract text content from a Markdown file"""
    try:
        with open(md_path, 'r', encoding='utf-8') as f:
            text = f.read()
        
        # Clean the text to remove problematic characters
        text = text.replace('\x00', '')
        text = text.replace('\ufffd', '')
        # Remove other non-printable characters except newlines and tabs
        text = ''.join(char for char in text if char.isprintable() or char in '\n\t')
        return text
    except Exception as e:
        raise RuntimeError(f"Error extracting text from {md_path}: {e}")


def extract_text_from_pdf(pdf_path):
    """Extract text content from a PDF file using PyMuPDF"""
    try:
        doc = fitz.open(pdf_path)
        text = ""
        for page_num, page in enumerate(doc):
            try:
                page_text = page.get_text()
                if page_text:
                    # Clean the text to remove problematic characters
                    page_text = page_text.replace('\x00', '')
                    page_text = page_text.replace('\ufffd', '')
                    # Remove other non-printable characters except newlines and tabs
                    page_text = ''.join(char for char in page_text if char.isprintable() or char in '\n\t')
                    text += page_text + "\n"
            except Exception as e:
                print(f"Warning: Could not extract text from page {page_num + 1}: {e}")
                continue
        doc.close()
        return text
    except Exception as e:
        raise RuntimeError(f"Error extracting text from {pdf_path}: {e}")


def chunk_text(text, chunk_size=None, overlap=None):
    """Split text into overlapping chunks for vector search"""
    if chunk_size is None:
        chunk_size = MAX_CHUNK_SIZE
    if overlap is None:
        overlap = CHUNK_OVERLAP
        
    if len(text) <= chunk_size:
        return [text]
    
    chunks = []
    start = 0
    
    while start < len(text):
        end = start + chunk_size
        
        # Try to break at sentence boundary
        if end < len(text):
            # Look for sentence endings near the chunk boundary
            for i in range(min(100, chunk_size // 4)):  # Look back up to 100 chars
                if end - i > start and text[end - i] in '.!?':
                    end = end - i + 1
                    break
        
        chunk = text[start:end].strip()
        if chunk:
            # Additional cleaning for each chunk
            chunk = chunk.replace('\x00', '').replace('\ufffd', '')
            # Ensure chunk only contains valid text
            chunk = ''.join(char for char in chunk if char.isprintable() or char in '\n\t ')
            # Replace multiple whitespaces with single space
            chunk = ' '.join(chunk.split())
            if chunk:  # Check again after cleaning
                chunks.append(chunk)
        
        start = end - overlap
        if start >= len(text):
            break
    
    return chunks


def get_file_hash(file_path):
    """Generate hash of file content for cache key"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


class PDFVectorStore:
    """FAISS-based vector store for PDF content"""
    
    def __init__(self, use_gpu=False, model_name='all-MiniLM-L6-v2', chunk_size=None, chunk_overlap=None):
        self.index = None
        self.chunks = []
        self.chunk_metadata = []  # Store filename and chunk info
        # Option to use faster model: 'paraphrase-MiniLM-L3-v2' (3x faster)
        self.embedder = SentenceTransformer(model_name)
        self.use_gpu = use_gpu
        self.embedding_dim = self.embedder.get_sentence_embedding_dimension()
        self.chunk_size = chunk_size if chunk_size is not None else MAX_CHUNK_SIZE
        self.chunk_overlap = chunk_overlap if chunk_overlap is not None else CHUNK_OVERLAP
        
        # Move embedder to GPU if available
        if use_gpu:
            try:
                self.embedder = self.embedder.cuda()
            except Exception as e:
                print(f"Warning: Could not move embedder to GPU: {e}")
                self.use_gpu = False
        
        # Initialize cache directory
        if ENABLE_EMBEDDING_CACHE:
            os.makedirs(CACHE_DIR, exist_ok=True)
    
    def add_documents(self, document_files):
        """Process PDF and Markdown files and add to vector store"""
        all_chunks = []
        all_metadata = []
        
        # Don't print this message here since spinner will show progress
        # print("Processing documents and extracting text...")
        
        # Process documents in parallel for better performance
        def process_single_document(doc_path):
            try:
                filename = os.path.basename(doc_path)
                
                # Check cache first
                if ENABLE_EMBEDDING_CACHE:
                    file_hash = get_file_hash(doc_path)
                    cache_key = f"{filename}_{file_hash}_{self.chunk_size}_{self.chunk_overlap}.pkl"
                    cache_path = os.path.join(CACHE_DIR, cache_key)
                    
                    if os.path.exists(cache_path):
                        try:
                            with open(cache_path, 'rb') as f:
                                cached_data = pickle.load(f)
                            print(f"  {filename}: Using cached embeddings ({cached_data['num_chunks']} chunks)")
                            return cached_data['filename'], cached_data['chunks'], cached_data['metadata'], cached_data['text_len'], cached_data['embeddings']
                        except Exception as e:
                            print(f"  Warning: Cache read failed for {filename}, regenerating: {e}")
                
                # Extract text based on file type
                if doc_path.lower().endswith('.pdf'):
                    text = extract_text_from_pdf(doc_path)
                elif doc_path.lower().endswith('.md'):
                    text = extract_text_from_markdown(doc_path)
                else:
                    raise ValueError(f"Unsupported file type: {doc_path}")
                
                chunks = chunk_text(text, self.chunk_size, self.chunk_overlap)
                
                metadata_list = []
                for i, chunk in enumerate(chunks):
                    metadata_list.append({
                        'filename': filename,
                        'chunk_id': i,
                        'total_chunks': len(chunks)
                    })
                
                return filename, chunks, metadata_list, len(text), None
            except Exception as e:
                print(f"Warning: Could not process {doc_path}: {e}")
                return None, None, None, None, None
        
        # Use ThreadPoolExecutor for parallel processing
        max_workers = min(len(document_files), 4)  # Limit concurrent threads
        cached_embeddings = []
        files_to_embed = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_doc = {executor.submit(process_single_document, doc_path): doc_path 
                           for doc_path in document_files}
            
            for future in as_completed(future_to_doc):
                filename, chunks, metadata_list, text_len, embeddings = future.result()
                if chunks:
                    if embeddings is not None:
                        # Use cached embeddings
                        cached_embeddings.extend(embeddings)
                    else:
                        # Need to generate embeddings
                        print(f"  {filename}: {len(chunks)} chunks from {text_len} characters")
                        files_to_embed.append((filename, chunks, metadata_list, text_len))
                    
                    all_chunks.extend(chunks)
                    all_metadata.extend(metadata_list)
        
        if not all_chunks:
            raise ValueError("No text content could be extracted from any PDF files")
        
        # Prepare embeddings
        all_embeddings = []
        
        # Add cached embeddings first
        if cached_embeddings:
            all_embeddings = cached_embeddings
            print(f"Loaded {len(cached_embeddings)} cached embeddings")
        
        # Generate embeddings for new files
        if files_to_embed:
            chunks_to_embed = []
            embed_metadata = []
            
            for filename, chunks, metadata_list, text_len in files_to_embed:
                chunks_to_embed.extend(chunks)
                embed_metadata.extend([(filename, i) for i in range(len(chunks))])
            
            print(f"Generating embeddings for {len(chunks_to_embed)} new chunks...")
            
            # Ensure all chunks are strings and handle any encoding issues
            clean_chunks = []
            for chunk in chunks_to_embed:
                if isinstance(chunk, str):
                    # Remove any problematic characters and ensure it's clean text
                    clean_chunk = chunk.strip().replace('\x00', '').replace('\ufffd', '')
                    if clean_chunk:  # Only add non-empty chunks
                        clean_chunks.append(clean_chunk)
                else:
                    clean_chunk = str(chunk).strip().replace('\x00', '').replace('\ufffd', '')
                    if clean_chunk:
                        clean_chunks.append(clean_chunk)
            
            if not clean_chunks:
                raise ValueError("No valid text chunks found after cleaning")
        
            try:
                # Try with batch processing to handle large datasets better
                embeddings = self.embedder.encode(
                    clean_chunks, 
                    show_progress_bar=True, 
                    convert_to_tensor=False,
                    batch_size=64,  # Increased batch size for faster processing
                    normalize_embeddings=True
                )
                
                # Cache the new embeddings
                if ENABLE_EMBEDDING_CACHE:
                    embed_idx = 0
                    for filename, chunks, metadata_list, text_len in files_to_embed:
                        file_embeddings = embeddings[embed_idx:embed_idx + len(chunks)]
                        embed_idx += len(chunks)
                        
                        # Find original file path
                        original_path = None
                        for doc_path in document_files:
                            if os.path.basename(doc_path) == filename:
                                original_path = doc_path
                                break
                        
                        if original_path:
                            file_hash = get_file_hash(original_path)
                            cache_key = f"{filename}_{file_hash}_{self.chunk_size}_{self.chunk_overlap}.pkl"
                            cache_path = os.path.join(CACHE_DIR, cache_key)
                            
                            cache_data = {
                                'filename': filename,
                                'chunks': chunks,
                                'metadata': metadata_list,
                                'text_len': text_len,
                                'embeddings': file_embeddings.tolist(),
                                'num_chunks': len(chunks)
                            }
                            
                            try:
                                with open(cache_path, 'wb') as f:
                                    pickle.dump(cache_data, f)
                            except Exception as e:
                                print(f"Warning: Could not cache embeddings for {filename}: {e}")
                
                all_embeddings.extend(embeddings.tolist())
                
            except Exception as e:
                print(f"Error during embedding generation: {e}")
                print("Trying alternative encoding method...")
                # Fallback: encode one by one to identify problematic chunks
                embeddings_list = []
                
                for i, chunk in enumerate(clean_chunks):
                    try:
                        # Ensure chunk is a proper string
                        if not isinstance(chunk, str) or not chunk.strip():
                            print(f"Skipping empty or invalid chunk {i}")
                            continue
                        
                        emb = self.embedder.encode([chunk], convert_to_tensor=False, normalize_embeddings=True)
                        embeddings_list.append(emb[0])
                    except Exception as chunk_error:
                        print(f"Skipping problematic chunk {i}: {chunk_error}")
                        if i < len(clean_chunks):
                            print(f"  Chunk preview: {clean_chunks[i][:50]}...")
                
                if not embeddings_list:
                    raise ValueError("Could not generate any embeddings")
                
                all_embeddings.extend(embeddings_list)
        
        # Convert all embeddings to numpy array
        embeddings = np.array(all_embeddings)
        
        # Create FAISS index
        dimension = embeddings.shape[1]
        num_vectors = embeddings.shape[0]
        
        # Use approximate search for large datasets
        if num_vectors > 1000:
            # Use IVF index for faster search on large datasets
            nlist = min(int(np.sqrt(num_vectors)), 100)  # number of clusters
            
            if self.use_gpu:
                # GPU index with IVF
                res = faiss.StandardGpuResources()
                quantizer = faiss.IndexFlatIP(dimension)
                index_cpu = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_INNER_PRODUCT)
                self.index = faiss.index_cpu_to_gpu(res, 0, index_cpu)
            else:
                # CPU index with IVF
                quantizer = faiss.IndexFlatIP(dimension)
                self.index = faiss.IndexIVFFlat(quantizer, dimension, nlist, faiss.METRIC_INNER_PRODUCT)
            
            # Train the index
            embeddings_np = embeddings.astype('float32')
            faiss.normalize_L2(embeddings_np)
            self.index.train(embeddings_np)
            self.index.add(embeddings_np)
            
            # Set search parameters for speed/accuracy tradeoff
            self.index.nprobe = min(10, nlist // 2)  # search this many clusters
        else:
            # Use exact search for small datasets
            if self.use_gpu:
                # GPU index
                res = faiss.StandardGpuResources()
                index_cpu = faiss.IndexFlatIP(dimension)
                self.index = faiss.index_cpu_to_gpu(res, 0, index_cpu)
            else:
                # CPU index
                self.index = faiss.IndexFlatIP(dimension)
            
            # Add embeddings to index
            embeddings_np = embeddings.astype('float32')
            faiss.normalize_L2(embeddings_np)  # Normalize for cosine similarity
            self.index.add(embeddings_np)
        
        self.chunks = all_chunks
        self.chunk_metadata = all_metadata
        
        print(f"Vector store ready with {len(all_chunks)} chunks")
    
    def search(self, query, k=TOP_K_CHUNKS):
        """Find most relevant chunks for a query"""
        if self.index is None:
            raise ValueError("No documents have been added to the vector store")
        
        clean_query = str(query).strip().replace('\x00', '').replace('\ufffd', '')
        # Additional cleaning for query
        clean_query = ''.join(char for char in clean_query if char.isprintable() or char in '\n\t ')
        clean_query = ' '.join(clean_query.split())
        
        if not clean_query:
            raise ValueError("Query is empty after cleaning")
        
        try:
            query_embedding = self.embedder.encode([clean_query], convert_to_tensor=False, normalize_embeddings=True)
            query_embedding = query_embedding.astype('float32')
            faiss.normalize_L2(query_embedding)
            
            # Ensure k doesn't exceed the number of chunks
            k = min(k, len(self.chunks))
            
            scores, indices = self.index.search(query_embedding, k)
            
            results = []
            for idx, i in enumerate(indices[0]):
                if i >= 0 and i < len(self.chunks):  # Valid index
                    results.append({
                        'chunk': self.chunks[i],
                        'score': float(scores[0][idx]),
                        'metadata': self.chunk_metadata[i]
                    })
            
            return results
        except Exception as e:
            print(f"Error during search: {e}")
            raise


def estimate_token_count(text):
    """
    Rough estimation of token count (approximately 4 characters per token for English text)
    """
    return len(text) // 4


def estimate_pdf_size_for_nova(pdf_paths):
    """
    Estimate if PDF files can fit in Nova model context window.
    Returns (can_fit, total_size_mb, estimated_tokens)
    """
    total_size_bytes = sum(os.path.getsize(path) for path in pdf_paths)
    total_size_mb = total_size_bytes / (1024 * 1024)
    
    # Nova models have document size limits:
    # - Individual documents: max 25MB
    # - Total documents: max 25MB when uploaded directly
    max_individual_size_mb = 25
    max_total_size_mb = 25
    
    # Check individual file sizes
    for pdf_path in pdf_paths:
        file_size_mb = os.path.getsize(pdf_path) / (1024 * 1024)
        if file_size_mb > max_individual_size_mb:
            return False, total_size_mb, 0
    
    # Check total size
    if total_size_mb > max_total_size_mb:
        return False, total_size_mb, 0
    
    # Rough estimate: 1MB ≈ 250,000 characters ≈ 62,500 tokens
    estimated_tokens = int(total_size_mb * 62500)
    
    return True, total_size_mb, estimated_tokens


def can_use_direct_pdf_submission(document_file_paths, model_id, question):
    """
    Determine if we should use direct PDF submission to Nova instead of FAISS.
    Returns True if:
    1. All files are PDFs (no markdown)
    2. Model supports document inputs
    3. Total size fits Nova limits
    4. Estimated tokens fit in model context window
    """
    # Only works with PDF files
    pdf_files = [path for path in document_file_paths if path.lower().endswith('.pdf')]
    if len(pdf_files) != len(document_file_paths):
        return False  # Has non-PDF files
    
    # Only works with models that support document inputs
    if not supports_document_input(model_id):
        # If using default model, check if we would auto-select a document-capable model
        if model_id == DEFAULT_BEDROCK_MODEL_ID:
            # Estimate if we would select a document-capable Nova model based on PDF size
            can_fit, total_size_mb, estimated_tokens = estimate_pdf_size_for_nova(pdf_files)
            if can_fit:
                # Would likely auto-select a document-capable Nova model
                return True
        return False
    
    # Check if PDFs fit Nova size limits
    can_fit, total_size_mb, estimated_tokens = estimate_pdf_size_for_nova(pdf_files)
    if not can_fit:
        return False
    
    # Estimate total prompt tokens (PDFs + question + overhead)
    question_tokens = estimate_token_count(question)
    overhead_tokens = 1000  # System prompt, formatting, etc.
    total_estimated_tokens = estimated_tokens + question_tokens + overhead_tokens
    
    # Check if it fits in the model's context window
    model_limit = MODEL_LIMITS.get(model_id, 128000)  # Default to micro if unknown
    
    return total_estimated_tokens <= (model_limit * 0.8)  # Use 80% of limit for safety


def submit_pdfs_directly_to_nova(bedrock_client, model_id, pdf_paths, question, inference_config):
    """
    Submit PDFs directly to Nova model using document support.
    Returns the response from the model.
    """
    # Prepare messages with PDF documents
    content = []
    
    # Add each PDF as a document
    for pdf_path in pdf_paths:
        filename = os.path.basename(pdf_path)
        sanitized_name = sanitize_document_name(filename)
        
        with open(pdf_path, 'rb') as f:
            pdf_bytes = f.read()
        
        content.append({
            "document": {
                "format": "pdf",
                "name": sanitized_name,
                "source": {
                    "bytes": pdf_bytes
                }
            }
        })
    
    # Add the question as text
    content.append({
        "text": question
    })
    
    messages = [{"role": "user", "content": content}]
    
    response = bedrock_client.converse(
        modelId=model_id,
        messages=messages,
        inferenceConfig=inference_config
    )
    
    return response


def select_model_for_context(context_text, preferred_model="us.amazon.nova-micro-v1:0"):
    """
    Select the appropriate Nova model based on context size.
    Returns the model ID that can handle the context size.
    """
    estimated_tokens = estimate_token_count(context_text)
    
    # Add some buffer for the response and system overhead
    required_tokens = estimated_tokens + 3000
    
    # Try models in order of preference (smallest to largest)
    model_order = [
        "us.amazon.nova-micro-v1:0",
        "us.amazon.nova-lite-v1:0", 
        "us.amazon.nova-premier-v1:0"
    ]
    
    for model_id in model_order:
        if required_tokens <= MODEL_LIMITS[model_id]:
            return model_id
    
    # If even the largest model can't handle it, return the largest and let it fail gracefully
    return "us.amazon.nova-premier-v1:0"


def is_nova_model(model_id):
    """Check if the model is a Nova model"""
    return "nova" in model_id.lower()


def supports_document_input(model_id):
    """Check if the model supports document inputs"""
    # Based on AWS documentation, only certain Nova models support document understanding
    # Nova Micro does not support document inputs
    document_supported_models = [
        "us.amazon.nova-lite-v1:0",
        "us.amazon.nova-pro-v1:0", 
        "us.amazon.nova-premier-v1:0"
    ]
    return model_id in document_supported_models


def estimate_max_chunks_for_model(model_id, base_prompt_size, avg_chunk_size):
    """
    Estimate maximum number of chunks that can fit in a model's context window.
    Returns the maximum number of chunks, or None if model limit is unknown.
    """
    if model_id in MODEL_LIMITS:
        max_tokens = MODEL_LIMITS[model_id]
    else:
        # For non-Nova models, assume a conservative limit
        max_tokens = 100000  # Conservative estimate for most models
    
    # Reserve tokens for response and overhead
    available_tokens = max_tokens - 3000
    
    # Convert to characters (rough estimate: 4 chars per token)
    available_chars = available_tokens * 4
    
    # Account for base prompt size
    remaining_chars = available_chars - base_prompt_size
    
    if remaining_chars <= 0:
        return 0
    
    # Calculate how many chunks can fit
    max_chunks = max(1, remaining_chars // avg_chunk_size)
    return max_chunks


def sanitize_document_name(filename):
    """
    Sanitizes a filename to be a valid document name for Bedrock.
    Keeps alphanumeric, hyphens, parentheses, square brackets. Replaces others with underscore.
    Limits length to 60 characters.
    """
    base_name = os.path.splitext(filename)[0]
    sanitized = re.sub(r'[^a-zA-Z0-9\-\(\)\[\]_]', '_', base_name)
    return sanitized[:60]


def get_document_files_details(input_path, recursive=False):
    """
    Collects PDF and Markdown file paths from a given file or directory.
    For directories, searches recursively if recursive=True, otherwise only top level.
    Returns a list of full file paths and the total size in bytes.
    """
    document_file_paths = []
    total_size_bytes = 0

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Error: Input path '{input_path}' not found.")

    if os.path.isfile(input_path):
        if input_path.lower().endswith((".pdf", ".md")):
            document_file_paths.append(input_path)
            total_size_bytes = os.path.getsize(input_path)
        else:
            raise ValueError(f"Error: Specified file '{input_path}' is not a PDF or Markdown file.")
    elif os.path.isdir(input_path):
        if recursive:
            # Search recursively through all subdirectories
            for root, dirs, files in os.walk(input_path):
                for file in files:
                    if file.lower().endswith((".pdf", ".md")):
                        file_path = os.path.join(root, file)
                        document_file_paths.append(file_path)
                        total_size_bytes += os.path.getsize(file_path)
        else:
            # Only look at files in the directory itself, not subdirectories
            for file in os.listdir(input_path):
                file_path = os.path.join(input_path, file)
                if os.path.isfile(file_path) and file.lower().endswith((".pdf", ".md")):
                    document_file_paths.append(file_path)
                    total_size_bytes += os.path.getsize(file_path)
    else:
        raise ValueError(f"Error: Input path '{input_path}' is not a valid file or directory.")

    if not document_file_paths:
        search_type = "recursively" if recursive else "in top-level directory"
        raise FileNotFoundError(f"Error: No PDF or Markdown files found {search_type} at '{input_path}'.")

    return document_file_paths, total_size_bytes




def main():
    parser = argparse.ArgumentParser(
        description="Ask questions about PDF documents using Amazon Bedrock (Nova Pro model)."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to a PDF/Markdown file or a folder containing PDF and Markdown files."
    )
    parser.add_argument(
        "question",
        type=str,
        help="The question to ask about the document(s)."
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=None,
        help="Number of most relevant chunks to retrieve for answering. Default: all chunks"
    )
    parser.add_argument(
        "--region",
        type=str,
        default=None,
        help=f"AWS region for Bedrock and S3. Default: region from AWS profile or {DEFAULT_AWS_REGION}"
    )
    parser.add_argument(
        "--profile",
        type=str,
        help="AWS profile name to use for credentials"
    )
    parser.add_argument(
        "--output-tokens",
        type=int,
        default=DEFAULT_OUTPUT_TOKENS,
        help=f"Maximum tokens for the model's output response. Default: {DEFAULT_OUTPUT_TOKENS}"
    )
    parser.add_argument(
        "--input-tokens",
        type=int,
        default=DEFAULT_INPUT_TOKENS,
        help=f"Informational: Target/assumed input token capacity. Not directly enforced by the API for document inputs. Default: {DEFAULT_INPUT_TOKENS}"
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=DEFAULT_TEMPERATURE,
        help=f"Temperature for model's response generation. Default: {DEFAULT_TEMPERATURE}"
    )
    parser.add_argument(
        "--top-p",
        type=float,
        default=DEFAULT_TOP_P,
        help=f"Top P for model's response generation. Default: {DEFAULT_TOP_P}"
    )
    parser.add_argument(
        "--model-id",
        type=str,
        default=DEFAULT_BEDROCK_MODEL_ID,
        help=f"The Bedrock model ID to use. Default: {DEFAULT_BEDROCK_MODEL_ID}"
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=MAX_CHUNK_SIZE,
        help=f"Size of text chunks for vector search. Larger = fewer chunks but less granular. Default: {MAX_CHUNK_SIZE}"
    )
    parser.add_argument(
        "--fast-embeddings",
        action="store_true",
        help="Use faster but potentially less accurate embedding model (3x speedup)"
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable embedding cache"
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search subdirectories recursively for PDF and Markdown files"
    )
    parser.add_argument(
        "--no-faiss",
        action="store_true",
        help="Disable FAISS vector store and use direct PDF submission when possible"
    )

    args = parser.parse_args()
    current_model_id = args.model_id
    response = None  # Initialize response variable

    # Start spinner for document discovery
    search_mode = "recursively" if args.recursive else "in top-level directory"
    spinner = Spinner(f"Scanning for documents {search_mode}")
    spinner.start()

    try:
        document_file_paths, total_size_bytes = get_document_files_details(args.path, recursive=args.recursive)
        spinner.stop()
    except (FileNotFoundError, ValueError) as e:
        spinner.stop()
        print(e, file=sys.stderr)
        sys.exit(1)

    num_documents = len(document_file_paths)
    pdf_count = sum(1 for path in document_file_paths if path.lower().endswith('.pdf'))
    md_count = sum(1 for path in document_file_paths if path.lower().endswith('.md'))
    print(f"Found {num_documents} document(s) ({pdf_count} PDF, {md_count} Markdown) {search_mode}, total size: {total_size_bytes / (1024*1024):.2f} MB.")

    # Initialize AWS session with optional profile
    session = boto3.Session(profile_name=args.profile) if args.profile else boto3.Session()
    
    # Determine region: command line arg > profile region > default
    region = args.region
    if not region:
        try:
            # Try to get region from session/profile
            region = session.region_name
            if region:
                print(f"Using region from AWS profile: {region}")
            else:
                region = DEFAULT_AWS_REGION
                print(f"No region found in profile, using default: {region}")
        except Exception:
            region = DEFAULT_AWS_REGION
            print(f"Could not determine region from profile, using default: {region}")
    else:
        print(f"Using region from command line: {region}")
    
    bedrock_client = session.client("bedrock-runtime", region_name=region)

    # Calculate chunk size and overlap
    chunk_size = args.chunk_size if args.chunk_size else MAX_CHUNK_SIZE
    chunk_overlap = max(200, chunk_size // 6)  # Keep proportional
    
    # Disable cache if requested
    if args.no_cache:
        global ENABLE_EMBEDDING_CACHE
        ENABLE_EMBEDDING_CACHE = False
    
    # By default use FAISS, only use direct PDF submission if --no-faiss is specified and criteria are met
    use_direct_pdf = False
    
    if args.no_faiss:
        # User wants to disable FAISS, check if direct PDF submission is possible
        use_direct_pdf = can_use_direct_pdf_submission(document_file_paths, current_model_id, args.question)
        if use_direct_pdf:
            print("FAISS disabled (--no-faiss specified), using direct PDF submission")
        else:
            print("FAISS disabled (--no-faiss specified) but direct PDF submission not possible, falling back to FAISS")
            use_direct_pdf = False
    else:
        print("Using FAISS vector store (default behavior)")
    
    # Initialize FAISS vector store only if needed
    vector_store = None
    if not use_direct_pdf:
        use_gpu = detect_gpu()
        embedding_model = 'paraphrase-MiniLM-L3-v2' if args.fast_embeddings else 'all-MiniLM-L6-v2'
        if args.fast_embeddings:
            print(f"Using fast embedding model: {embedding_model}")
        
        # Print FAISS detection message only when actually using FAISS
        if not use_gpu:
            print("FAISS CPU version detected, using CPU-only processing")
        
        vector_store = PDFVectorStore(
            use_gpu=use_gpu, 
            model_name=embedding_model,
            chunk_size=chunk_size,
            chunk_overlap=chunk_overlap
        )
    
    if use_direct_pdf:
        # Use direct PDF submission to Nova
        inference_config = {
            "maxTokens": args.output_tokens,
            "temperature": args.temperature,
            "topP": args.top_p
        }
        
        pdf_files = [path for path in document_file_paths if path.lower().endswith('.pdf')]
        can_fit, total_size_mb, estimated_tokens = estimate_pdf_size_for_nova(pdf_files)
        
        # Auto-select appropriate Nova model if using default
        if args.model_id == DEFAULT_BEDROCK_MODEL_ID:
            # Select model based on estimated tokens, but only document-capable models
            question_tokens = estimate_token_count(args.question)
            overhead_tokens = 1000
            total_estimated_tokens = estimated_tokens + question_tokens + overhead_tokens
            
            # Skip Nova Micro since it doesn't support documents
            if total_estimated_tokens <= MODEL_LIMITS["us.amazon.nova-lite-v1:0"] * 0.8:
                current_model_id = "us.amazon.nova-lite-v1:0"
            elif total_estimated_tokens <= MODEL_LIMITS["us.amazon.nova-pro-v1:0"] * 0.8:
                current_model_id = "us.amazon.nova-pro-v1:0"
            else:
                current_model_id = "us.amazon.nova-premier-v1:0"
            
            print(f"Auto-selected {current_model_id} for direct PDF submission")
        
        print(f"Direct PDF submission: {len(pdf_files)} PDF(s), {total_size_mb:.2f} MB, ~{estimated_tokens} tokens")
        print(f"Using direct PDF submission to Nova model (bypassing FAISS)")
        
        try:
            # Start spinner for PDF submission
            spinner = Spinner(f"Submitting {len(pdf_files)} PDF(s) directly to {current_model_id}")
            spinner.start()
            
            response = submit_pdfs_directly_to_nova(
                bedrock_client, current_model_id, pdf_files, args.question, inference_config
            )
            
            spinner.stop()
            print()  # Add a blank line after spinner stops
        except (ClientError, ReadTimeoutError) as e:
            error_code = e.response.get("Error", {}).get("Code")
            error_message = e.response.get("Error", {}).get("Message", str(e))
            
            # Stop spinner on error
            if 'spinner' in locals():
                spinner.stop()
            
            # Handle different types of errors
            if isinstance(e, ReadTimeoutError):
                print(f"Direct PDF submission timed out, falling back to FAISS approach...")
                use_direct_pdf = False
            elif isinstance(e, ClientError):
                error_code = e.response.get("Error", {}).get("Code")
                error_message = e.response.get("Error", {}).get("Message", str(e))
                
                # If direct submission fails, fall back to FAISS
                if (error_code == "ValidationException" and 
                    ("Input is too long" in error_message or 
                     "too long for requested model" in error_message or
                     "doesn't support documents" in error_message)):
                    print(f"Direct PDF submission failed ({error_message}), falling back to FAISS approach...")
                    use_direct_pdf = False
                else:
                    raise e
            else:
                raise e
    
    if not use_direct_pdf:
        # Use FAISS-based approach
        try:
            # Build vector store from documents with spinner
            spinner = Spinner("Processing documents and building vector store")
            spinner.start()
            
            vector_store.add_documents(document_file_paths)
            
            spinner.stop()
            
            # Search for relevant chunks
            print(f"\nSearching for relevant content for: '{args.question}'")
            # Use all chunks if top_k is None, otherwise use specified number
            search_k = args.top_k if args.top_k is not None else len(vector_store.chunks)
            relevant_chunks = vector_store.search(args.question, k=search_k)
            
            if not relevant_chunks:
                print("No relevant content found in the documents.", file=sys.stderr)
                sys.exit(1)
            
            # Display found chunks
            print(f"\nFound {len(relevant_chunks)} relevant chunks:")
            for i, result in enumerate(relevant_chunks):
                metadata = result['metadata']
                print(f"  {i+1}. {metadata['filename']} (chunk {metadata['chunk_id']+1}/{metadata['total_chunks']}) - Score: {result['score']:.3f}")
            
            # Combine relevant chunks into context
            context_parts = []
            for result in relevant_chunks:
                metadata = result['metadata']
                chunk_header = f"[From {metadata['filename']}, section {metadata['chunk_id']+1}]"
                context_parts.append(f"{chunk_header}\n{result['chunk']}")
            
            combined_context = "\n\n---\n\n".join(context_parts)
            
            # Create prompt with context and question
            prompt = f"""Based on the following document excerpts, please answer the question.

Document excerpts:
{combined_context}

Question: {args.question}

Please provide a comprehensive answer based on the information in the document excerpts above."""

            # Auto-select model based on context size only if using default model
            if args.model_id == DEFAULT_BEDROCK_MODEL_ID:
                selected_model = select_model_for_context(prompt, current_model_id)
                if selected_model != current_model_id:
                    print(f"\nAuto-selecting model {selected_model} based on context size ({len(prompt)} characters)")
                    current_model_id = selected_model
            else:
                print(f"\nUsing user-specified model: {current_model_id}")

            messages = [{"role": "user", "content": [{"text": prompt}]}]

            inference_config = {
                "maxTokens": args.output_tokens,
                "temperature": args.temperature,
                "topP": args.top_p
            }

            print(f"\nSending request to Amazon Bedrock ({current_model_id})...")
            print(f"Context size: {len(combined_context)} characters")
            print(f"Total prompt size: {len(prompt)} characters (~{estimate_token_count(prompt)} tokens)")
            
            response = None
            
            try:
                response = bedrock_client.converse(
                    modelId=current_model_id,
                    messages=messages,
                    inferenceConfig=inference_config
                )
                
            except (ClientError, ReadTimeoutError) as e:
                error_code = e.response.get("Error", {}).get("Code")
                error_message = e.response.get("Error", {}).get("Message", str(e))
                
                # Handle timeout errors
                if isinstance(e, ReadTimeoutError):
                    print(f"Request timed out with model {current_model_id}. This may indicate the document is too large or complex.")
                    print("Consider using --no-faiss with smaller documents or breaking large documents into smaller chunks.")
                    raise e
                
                # Check if it's an "input too long" error
                elif (isinstance(e, ClientError) and 
                      error_code == "ValidationException" and 
                      ("Input is too long" in error_message or "too long for requested model" in error_message)):
                
                    if is_nova_model(current_model_id):
                        # For Nova models, try fallback to larger Nova models
                        print(f"Model {current_model_id} cannot handle input size, trying larger Nova models...")
                        
                        fallback_models = [
                            "us.amazon.nova-lite-v1:0",
                            "us.amazon.nova-premier-v1:0"
                        ]
                        
                        # Remove current model if it's already in the list
                        fallback_models = [m for m in fallback_models if m != current_model_id]
                    
                        for attempt_model in fallback_models:
                            try:
                                print(f"\nFalling back to larger Nova model: {attempt_model}")
                                
                                response = bedrock_client.converse(
                                    modelId=attempt_model,
                                    messages=messages,
                                    inferenceConfig=inference_config
                                )
                                
                                current_model_id = attempt_model
                                break
                                
                            except (ClientError, ReadTimeoutError) as fallback_e:
                                if isinstance(fallback_e, ReadTimeoutError):
                                    print(f"Timeout with {attempt_model}, trying next model...")
                                    continue
                                elif isinstance(fallback_e, ClientError):
                                    fallback_error_code = fallback_e.response.get("Error", {}).get("Code")
                                    fallback_error_message = fallback_e.response.get("Error", {}).get("Message", str(fallback_e))
                                    
                                    if (fallback_error_code == "ValidationException" and 
                                        ("Input is too long" in fallback_error_message or "too long for requested model" in fallback_error_message)):
                                        continue
                                    else:
                                        raise fallback_e
                                else:
                                    raise fallback_e
                        
                        if response is None:
                            raise e  # All Nova models failed
                            
                    else:
                        # For non-Nova models, reduce top-k and retry
                        print(f"Non-Nova model {current_model_id} cannot handle input size, reducing chunks...")
                    
                        # Calculate average chunk size
                        if relevant_chunks:
                            avg_chunk_size = sum(len(chunk['chunk']) for chunk in relevant_chunks) // len(relevant_chunks)
                            base_prompt_template = f"""Based on the following document excerpts, please answer the question.

Document excerpts:
{{CONTEXT}}

Question: {args.question}

Please provide a comprehensive answer based on the information in the document excerpts above."""
                            base_prompt_size = len(base_prompt_template) - len("{CONTEXT}")
                            
                            # Estimate maximum chunks that can fit
                            max_chunks = estimate_max_chunks_for_model(current_model_id, base_prompt_size, avg_chunk_size)
                            
                            if max_chunks > 0 and max_chunks < len(relevant_chunks):
                                print(f"Reducing from {len(relevant_chunks)} to {max_chunks} chunks to fit context window")
                                
                                # Use only the top-scoring chunks
                                reduced_chunks = relevant_chunks[:max_chunks]
                                
                                # Rebuild context with reduced chunks
                                context_parts = []
                                for result in reduced_chunks:
                                    metadata = result['metadata']
                                    chunk_header = f"[From {metadata['filename']}, section {metadata['chunk_id']+1}]"
                                    context_parts.append(f"{chunk_header}\n{result['chunk']}")
                                
                                reduced_context = "\n\n---\n\n".join(context_parts)
                                
                                # Create new prompt with reduced context
                                reduced_prompt = f"""Based on the following document excerpts, please answer the question.

Document excerpts:
{reduced_context}

Question: {args.question}

Please provide a comprehensive answer based on the information in the document excerpts above."""

                                reduced_messages = [{"role": "user", "content": [{"text": reduced_prompt}]}]
                                
                                print(f"Retrying with reduced context: {len(reduced_context)} characters (~{estimate_token_count(reduced_prompt)} tokens)")
                                
                                try:
                                    response = bedrock_client.converse(
                                        modelId=current_model_id,
                                        messages=reduced_messages,
                                        inferenceConfig=inference_config
                                    )
                                except (ClientError, ReadTimeoutError) as retry_e:
                                    # If it still fails, raise the original error
                                    if isinstance(retry_e, ReadTimeoutError):
                                        print(f"Request still timed out after reducing context size.")
                                    raise e
                            else:
                                raise e  # Can't reduce further
                        else:
                            raise e  # No chunks to work with
                else:
                    # For other errors, don't retry
                    raise e

        except Exception as e:
            print(f"\nAn unexpected error occurred during FAISS processing: {e}", file=sys.stderr)
            sys.exit(1)

    # Extract and display response (works for both direct PDF and FAISS approaches)
    if response:
        response_text = response['output']['message']['content'][0]['text']
        print("\n[Model Response]")
        print(response_text)
    else:
        print("No response received from the model.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
