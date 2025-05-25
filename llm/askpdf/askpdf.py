#!/usr/bin/env python3

import argparse
import boto3
import os
import json
import uuid
import re
import sys
import tempfile
import shutil
from botocore.exceptions import ClientError
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# FAISS and embedding imports
try:
    import faiss
except ImportError:
    print(f"Error: FAISS not installed. Run: pip install faiss-gpu sentence-transformers pymupdf", file=sys.stderr)
    print(f"For CPU-only FAISS: pip install faiss-cpu sentence-transformers pymupdf", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error: FAISS installation appears corrupted. Try reinstalling:", file=sys.stderr)
    print(f"  pip uninstall faiss-gpu faiss-cpu", file=sys.stderr)
    print(f"  pip install faiss-cpu sentence-transformers pymupdf", file=sys.stderr)
    print(f"Original error: {e}", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
    from sentence_transformers import SentenceTransformer
    try:
        import fitz  # PyMuPDF
        USE_PYMUPDF = True
    except ImportError:
        try:
            from PyPDF2 import PdfReader
            USE_PYMUPDF = False
            print("Warning: PyMuPDF not installed, falling back to slower PyPDF2. For better performance, run: pip install pymupdf")
        except ImportError:
            print(f"Error: No PDF library installed. Run: pip install pymupdf sentence-transformers", file=sys.stderr)
            sys.exit(1)
except ImportError as e:
    print(f"Error: Required packages not installed. Run: pip install pymupdf sentence-transformers", file=sys.stderr)
    sys.exit(1)

# --- Configuration Constants ---
# Default Model ID for Amazon Titan
DEFAULT_BEDROCK_MODEL_ID = "us.amazon.nova-pro-v1:0"
# AWS Region for Bedrock and S3 services
AWS_REGION = "us-east-1" # You can change this if needed

# Document processing limits - now used for chunking strategy
MAX_CHUNK_SIZE = 1000  # Characters per chunk for vector search
CHUNK_OVERLAP = 200    # Overlap between chunks
TOP_K_CHUNKS = 5       # Number of most relevant chunks to send to LLM

# Bedrock inference parameters
DEFAULT_INPUT_TOKENS = 30000  # Placeholder for typical model context window
DEFAULT_OUTPUT_TOKENS = 2048
DEFAULT_TOP_P = 0.9
DEFAULT_TEMPERATURE = 0.2


def detect_gpu():
    """Detect if NVIDIA GPU is available for FAISS"""
    try:
        import torch
        if torch.cuda.is_available():
            gpu_count = torch.cuda.device_count()
            gpu_name = torch.cuda.get_device_name(0)
            print(f"GPU detected: {gpu_name} (using GPU-accelerated FAISS)")
            return True
    except ImportError:
        pass
    
    # Alternative check using nvidia-ml-py
    try:
        import pynvml
        pynvml.nvmlInit()
        gpu_count = pynvml.nvmlDeviceGetCount()
        if gpu_count > 0:
            print(f"GPU detected: {gpu_count} NVIDIA GPU(s) (using GPU-accelerated FAISS)")
            return True
    except ImportError:
        pass
    
    print("No GPU detected, using CPU FAISS")
    return False


def extract_text_from_pdf(pdf_path):
    """Extract text content from a PDF file"""
    if USE_PYMUPDF:
        # Use PyMuPDF (much faster)
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
    else:
        # Fallback to PyPDF2
        try:
            reader = PdfReader(pdf_path)
            text = ""
            for page_num, page in enumerate(reader.pages):
                try:
                    page_text = page.extract_text()
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
            return text
        except Exception as e:
            raise RuntimeError(f"Error extracting text from {pdf_path}: {e}")


def chunk_text(text, chunk_size=MAX_CHUNK_SIZE, overlap=CHUNK_OVERLAP):
    """Split text into overlapping chunks for vector search"""
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


class PDFVectorStore:
    """FAISS-based vector store for PDF content"""
    
    def __init__(self, use_gpu=False):
        self.index = None
        self.chunks = []
        self.chunk_metadata = []  # Store filename and chunk info
        self.embedder = SentenceTransformer('all-MiniLM-L6-v2')
        self.use_gpu = use_gpu
        
        # Move embedder to GPU if available
        if use_gpu:
            try:
                self.embedder = self.embedder.cuda()
            except Exception as e:
                print(f"Warning: Could not move embedder to GPU: {e}")
                self.use_gpu = False
    
    def add_documents(self, pdf_files):
        """Process PDF files and add to vector store"""
        all_chunks = []
        all_metadata = []
        
        print("Processing PDFs and extracting text...")
        
        # Process PDFs in parallel for better performance
        def process_single_pdf(pdf_path):
            try:
                text = extract_text_from_pdf(pdf_path)
                chunks = chunk_text(text)
                filename = os.path.basename(pdf_path)
                
                metadata_list = []
                for i, chunk in enumerate(chunks):
                    metadata_list.append({
                        'filename': filename,
                        'chunk_id': i,
                        'total_chunks': len(chunks)
                    })
                
                return filename, chunks, metadata_list, len(text)
            except Exception as e:
                print(f"Warning: Could not process {pdf_path}: {e}")
                return None, None, None, None
        
        # Use ThreadPoolExecutor for parallel processing
        max_workers = min(len(pdf_files), 4)  # Limit concurrent threads
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_pdf = {executor.submit(process_single_pdf, pdf_path): pdf_path 
                           for pdf_path in pdf_files}
            
            for future in as_completed(future_to_pdf):
                filename, chunks, metadata_list, text_len = future.result()
                if chunks:
                    print(f"  {filename}: {len(chunks)} chunks from {text_len} characters")
                    all_chunks.extend(chunks)
                    all_metadata.extend(metadata_list)
        
        if not all_chunks:
            raise ValueError("No text content could be extracted from any PDF files")
        
        print(f"Generating embeddings for {len(all_chunks)} chunks...")
        # Ensure all chunks are strings and handle any encoding issues
        clean_chunks = []
        for chunk in all_chunks:
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
                batch_size=32,
                normalize_embeddings=True
            )
        except Exception as e:
            print(f"Error during embedding generation: {e}")
            print("Trying alternative encoding method...")
            # Fallback: encode one by one to identify problematic chunks
            embeddings_list = []
            valid_chunks = []
            valid_metadata = []
            
            for i, chunk in enumerate(clean_chunks):
                try:
                    # Ensure chunk is a proper string
                    if not isinstance(chunk, str) or not chunk.strip():
                        print(f"Skipping empty or invalid chunk {i}")
                        continue
                    
                    emb = self.embedder.encode([chunk], convert_to_tensor=False, normalize_embeddings=True)
                    embeddings_list.append(emb[0])
                    valid_chunks.append(chunk)
                    valid_metadata.append(all_metadata[i])
                except Exception as chunk_error:
                    print(f"Skipping problematic chunk {i}: {chunk_error}")
                    if i < len(clean_chunks):
                        print(f"  Chunk preview: {clean_chunks[i][:50]}...")
            
            if not embeddings_list:
                raise ValueError("Could not generate any embeddings")
            
            embeddings = np.array(embeddings_list)
            # Update chunks and metadata to only include valid ones
            all_chunks = valid_chunks
            all_metadata = valid_metadata
        
        # Create FAISS index
        dimension = embeddings.shape[1]
        
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


def sanitize_document_name(filename):
    """
    Sanitizes a filename to be a valid document name for Bedrock.
    Keeps alphanumeric, hyphens, parentheses, square brackets. Replaces others with underscore.
    Limits length to 60 characters.
    """
    base_name = os.path.splitext(filename)[0]
    sanitized = re.sub(r'[^a-zA-Z0-9\-\(\)\[\]_]', '_', base_name)
    return sanitized[:60]


def get_pdf_files_details(input_path):
    """
    Collects PDF file paths from a given file or directory.
    Returns a list of full file paths and the total size in bytes.
    """
    pdf_file_paths = []
    total_size_bytes = 0

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Error: Input path '{input_path}' not found.")

    if os.path.isfile(input_path):
        if input_path.lower().endswith(".pdf"):
            pdf_file_paths.append(input_path)
            total_size_bytes = os.path.getsize(input_path)
        else:
            raise ValueError(f"Error: Specified file '{input_path}' is not a PDF.")
    elif os.path.isdir(input_path):
        for root, _, files in os.walk(input_path):
            for file in files:
                if file.lower().endswith(".pdf"):
                    full_path = os.path.join(root, file)
                    pdf_file_paths.append(full_path)
                    total_size_bytes += os.path.getsize(full_path)
    else:
        raise ValueError(f"Error: Input path '{input_path}' is not a valid file or directory.")

    if not pdf_file_paths:
        raise FileNotFoundError(f"Error: No PDF files found at '{input_path}'.")

    return pdf_file_paths, total_size_bytes




def main():
    parser = argparse.ArgumentParser(
        description="Ask questions about PDF documents using Amazon Bedrock (Nova Pro model)."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to a PDF file or a folder containing PDF files."
    )
    parser.add_argument(
        "question",
        type=str,
        help="The question to ask about the PDF(s)."
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=TOP_K_CHUNKS,
        help=f"Number of most relevant chunks to retrieve for answering. Default: {TOP_K_CHUNKS}"
    )
    parser.add_argument(
        "--region",
        type=str,
        default=AWS_REGION,
        help=f"AWS region for Bedrock and S3. Default: {AWS_REGION}"
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

    args = parser.parse_args()
    current_model_id = args.model_id

    try:
        pdf_file_paths, total_size_bytes = get_pdf_files_details(args.path)
    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    num_documents = len(pdf_file_paths)
    print(f"Found {num_documents} PDF document(s), total size: {total_size_bytes / (1024*1024):.2f} MB.")

    # Initialize AWS session with optional profile
    session = boto3.Session(profile_name=args.profile) if args.profile else boto3.Session()
    bedrock_client = session.client("bedrock-runtime", region_name=args.region)

    # Initialize FAISS vector store
    use_gpu = detect_gpu()
    vector_store = PDFVectorStore(use_gpu=use_gpu)

    try:
        # Build vector store from PDFs
        vector_store.add_documents(pdf_file_paths)
        
        # Search for relevant chunks
        print(f"\nSearching for relevant content for: '{args.question}'")
        relevant_chunks = vector_store.search(args.question, k=args.top_k)
        
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

        messages = [{"role": "user", "content": [{"text": prompt}]}]

        inference_config = {
            "maxTokens": args.output_tokens,
            "temperature": args.temperature,
            "topP": args.top_p
        }

        print(f"\nSending request to Amazon Bedrock ({current_model_id})...")
        print(f"Context size: {len(combined_context)} characters")
        
        response = bedrock_client.converse(
            modelId=current_model_id,
            messages=messages,
            inferenceConfig=inference_config
        )

        response_text = response['output']['message']['content'][0]['text']
        print("\n[Model Response]")
        print(response_text)

    except ClientError as e:
        error_code = e.response.get("Error", {}).get("Code")
        error_message = e.response.get("Error", {}).get("Message", str(e))
        
        base_error_text = f"A Bedrock API error occurred with model {current_model_id}"
        specific_error_info = f"Bedrock Error Code: {error_code}, Message: {error_message}"

        full_message = f"\n{base_error_text}.\n{specific_error_info}"

        if error_code == "ValidationException" and "Input is too long" in error_message:
            full_message += (
                "\n\n[Additional Diagnostics]:"
                f"\nThe context sent to the model exceeded the token processing limit."
                f"\nModel used: {current_model_id}"
                "\n\nPossible solutions:"
                "\n1. Reduce --top-k to retrieve fewer chunks (currently using top-{args.top_k})"
                "\n2. Try a different Nova model with a larger context window:"
                "\n   --model-id us.amazon.nova-premier-v1:0  (largest context window)"
                "\n   --model-id us.amazon.nova-lite-v1:0     (smaller but more efficient)"
                "\n3. The retrieved chunks may contain dense information"
            )
        
        print(full_message, file=sys.stderr)
        sys.exit(1)
    except Exception as e: # General fallback for non-ClientError exceptions
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # No cleanup needed for FAISS - everything is in memory
        pass

if __name__ == "__main__":
    main()
