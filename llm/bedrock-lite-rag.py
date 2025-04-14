"""
Enhanced Serverless RAG Application using AWS Bedrock and S3
===========================================================

This implementation provides a complete solution for:
1. Setting up the required AWS resources
2. Processing PDFs from S3 (including large volumes)
3. Creating temporary vector stores with checkpointing
4. Distributed embedding generation for large documents
5. Memory-optimized querying with FAISS
6. Tearing down resources when done

Prerequisites:
- AWS CLI configured
- boto3
- Required Python libraries
"""

import os
import boto3
import json
import uuid
import time
import logging
import gc
import math
from typing import List, Dict, Any, Optional, Tuple
import tempfile
import shutil
from datetime import datetime, timedelta
from concurrent.futures import ThreadPoolExecutor

# For PDF processing
import PyPDF2
from langchain.text_splitter import RecursiveCharacterTextSplitter

# For vector operations
import numpy as np
from faiss import IndexFlatL2, write_index, read_index
import pickle
import mmap

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# AWS Service clients
s3_client = boto3.client('s3')
bedrock_client = boto3.client('bedrock-runtime')
lambda_client = boto3.client('lambda')

# Configuration
DEFAULT_BUCKET_NAME = "serverless-rag-pdf-storage"
TEMP_DIR = "temp_vectors"
CHECKPOINT_DIR = "checkpoints"
SESSION_TIMEOUT = 60 * 60  # 1 hour in seconds
CLAUDE_MODEL_ID = "anthropic.claude-3-sonnet-20240229-v1:0"  # Update as needed
EMBEDDING_MODEL_ID = "amazon.titan-embed-text-v1"  # Update as needed
CHUNK_SIZE = 1000
CHUNK_OVERLAP = 200
BATCH_SIZE = 10  # Number of PDFs to process in a batch
EMBEDDING_BATCH_SIZE = 20  # Number of text chunks to embed in a batch
LARGE_INDEX_THRESHOLD = 100000  # When to use memory-mapped files

class EnhancedServerlessRAG:
    def __init__(self, bucket_name: str = DEFAULT_BUCKET_NAME, session_id: str = None):
        """Initialize the RAG application."""
        self.bucket_name = bucket_name
        self.session_id = session_id or str(uuid.uuid4())
        self.vector_store = None
        self.index_to_text_map = {}
        self.last_activity = time.time()
        self.processed_files = set()
        self.use_mmap = False
        
        # Create temporary directories for storing vector data and checkpoints
        os.makedirs(TEMP_DIR, exist_ok=True)
        os.makedirs(os.path.join(TEMP_DIR, CHECKPOINT_DIR), exist_ok=True)
        
        # Check if bucket exists, if not create it
        self._ensure_bucket_exists()
    
    def _ensure_bucket_exists(self):
        """Create S3 bucket if it doesn't exist."""
        try:
            s3_client.head_bucket(Bucket=self.bucket_name)
            logger.info(f"Bucket {self.bucket_name} already exists")
        except:
            logger.info(f"Creating bucket {self.bucket_name}")
            region = boto3.session.Session().region_name
            
            # Different API call for us-east-1
            if region == 'us-east-1':
                s3_client.create_bucket(Bucket=self.bucket_name)
            else:
                s3_client.create_bucket(
                    Bucket=self.bucket_name,
                    CreateBucketConfiguration={'LocationConstraint': region}
                )
            
            # Enable versioning for data protection
            s3_client.put_bucket_versioning(
                Bucket=self.bucket_name,
                VersioningConfiguration={'Status': 'Enabled'}
            )
            
            # Create necessary folders
            s3_client.put_object(Bucket=self.bucket_name, Key=f"{CHECKPOINT_DIR}/")
            s3_client.put_object(Bucket=self.bucket_name, Key="temp/")
    
    def upload_pdf(self, file_path: str, key: Optional[str] = None) -> str:
        """Upload a PDF to S3.
        
        Args:
            file_path: Local path to the PDF file
            key: Optional S3 key, if None will use the filename
            
        Returns:
            S3 key of the uploaded file
        """
        if key is None:
            key = os.path.basename(file_path)
        
        logger.info(f"Uploading {file_path} to S3 bucket {self.bucket_name} as {key}")
        s3_client.upload_file(file_path, self.bucket_name, key)
        return key
    
    def list_pdfs(self) -> List[str]:
        """List all PDFs in the S3 bucket."""
        response = s3_client.list_objects_v2(Bucket=self.bucket_name)
        
        if 'Contents' not in response:
            return []
        
        return [obj['Key'] for obj in response['Contents'] 
                if obj['Key'].lower().endswith('.pdf')]
    
    def extract_text_from_pdf(self, s3_key: str) -> str:
        """Extract text from a PDF stored in S3."""
        with tempfile.NamedTemporaryFile(suffix='.pdf', delete=False) as temp_file:
            temp_path = temp_file.name
            
        try:
            logger.info(f"Downloading {s3_key} from S3")
            s3_client.download_file(self.bucket_name, s3_key, temp_path)
            
            logger.info(f"Extracting text from {s3_key}")
            with open(temp_path, 'rb') as f:
                pdf_reader = PyPDF2.PdfReader(f)
                text = ""
                
                # Process pages in batches to manage memory
                total_pages = len(pdf_reader.pages)
                batch_size = min(50, total_pages)  # Process 50 pages at a time
                
                for batch_start in range(0, total_pages, batch_size):
                    batch_end = min(batch_start + batch_size, total_pages)
                    logger.info(f"Processing pages {batch_start+1}-{batch_end} of {total_pages}")
                    
                    for page_num in range(batch_start, batch_end):
                        page = pdf_reader.pages[page_num]
                        text += page.extract_text() + "\n\n"
                    
                    # Force garbage collection between batches
                    gc.collect()
            
            return text
        except Exception as e:
            logger.error(f"Error extracting text from {s3_key}: {str(e)}")
            return f"Error processing document {s3_key}: {str(e)}"
        finally:
            # Clean up the temporary file
            if os.path.exists(temp_path):
                os.remove(temp_path)
    
    def chunk_text(self, text: str) -> List[str]:
        """Split text into manageable chunks for embedding."""
        text_splitter = RecursiveCharacterTextSplitter(
            chunk_size=CHUNK_SIZE,
            chunk_overlap=CHUNK_OVERLAP,
            length_function=len
        )
        return text_splitter.split_text(text)
    
    def get_embeddings(self, texts: List[str]) -> List[List[float]]:
        """Generate embeddings for text chunks using AWS Bedrock."""
        embeddings = []
        
        # Process in batches to manage memory and API rate limits
        for batch_idx in range(0, len(texts), EMBEDDING_BATCH_SIZE):
            batch_end = min(batch_idx + EMBEDDING_BATCH_SIZE, len(texts))
            batch = texts[batch_idx:batch_end]
            
            logger.info(f"Generating embeddings for chunks {batch_idx+1}-{batch_end} of {len(texts)}")
            
            batch_embeddings = []
            for text in batch:
                try:
                    # Format the request for the embedding model
                    request_body = json.dumps({
                        "inputText": text
                    })
                    
                    response = bedrock_client.invoke_model(
                        modelId=EMBEDDING_MODEL_ID,
                        contentType="application/json",
                        accept="application/json",
                        body=request_body
                    )
                    
                    response_body = json.loads(response['body'].read())
                    embedding_vector = response_body['embedding']
                    batch_embeddings.append(embedding_vector)
                    
                except Exception as e:
                    logger.error(f"Error generating embedding: {str(e)}")
                    # Return a zero vector as fallback
                    batch_embeddings.append([0.0] * 1536)  # Default size for many embedding models
            
            embeddings.extend(batch_embeddings)
            
            # Force garbage collection between batches
            gc.collect()
        
        return embeddings
    
    def get_embeddings_distributed(self, texts: List[str]) -> List[List[float]]:
        """Generate embeddings using distributed Lambda functions for large volumes."""
        if len(texts) <= EMBEDDING_BATCH_SIZE:
            # For small batches, use the regular method
            return self.get_embeddings(texts)
        
        # For larger batches, distribute the work
        embeddings = []
        
        # Split texts into chunks that can be processed by Lambda
        chunks = [texts[i:i+EMBEDDING_BATCH_SIZE] for i in range(0, len(texts), EMBEDDING_BATCH_SIZE)]
        logger.info(f"Distributing embedding generation across {len(chunks)} Lambda invocations")
        
        # Process each chunk with Lambda
        for i, chunk in enumerate(chunks):
            logger.info(f"Processing chunk {i+1}/{len(chunks)} with Lambda")
            
            try:
                # Invoke Lambda function for embedding generation
                lambda_payload = {
                    "texts": chunk,
                    "model_id": EMBEDDING_MODEL_ID
                }
                
                response = lambda_client.invoke(
                    FunctionName="bedrock-embedding-generator",
                    InvocationType='RequestResponse',
                    Payload=json.dumps(lambda_payload)
                )
                
                # Parse response
                response_payload = json.loads(response['Payload'].read())
                if 'embeddings' in response_payload:
                    embeddings.extend(response_payload['embeddings'])
                else:
                    # Handle error by adding zero vectors
                    logger.error(f"Lambda error: {response_payload.get('error', 'Unknown error')}")
                    embeddings.extend([[0.0] * 1536] * len(chunk))
            except Exception as e:
                logger.error(f"Error invoking Lambda: {str(e)}")
                # Add zero vectors as fallback
                embeddings.extend([[0.0] * 1536] * len(chunk))
        
        return embeddings
    
    def build_index(self, pdfs: Optional[List[str]] = None) -> None:
        """Process PDFs and build a FAISS index for vector search.
        
        Args:
            pdfs: List of S3 keys for PDFs to process. If None, process all PDFs.
        """
        # If no PDFs specified, get all PDFs from the bucket
        if pdfs is None:
            pdfs = self.list_pdfs()
        
        if not pdfs:
            logger.warning("No PDFs found to process")
            return
        
        # Check if we have a large number of PDFs
        if len(pdfs) > BATCH_SIZE:
            logger.info(f"Large PDF set detected ({len(pdfs)} files). Using batch processing.")
            self.process_large_pdfs(pdfs)
            return
        
        all_chunks = []
        chunk_to_source = {}
        chunk_idx = 0
        
        # Process each PDF
        for pdf_key in pdfs:
            if pdf_key in self.processed_files:
                logger.info(f"Skipping already processed file: {pdf_key}")
                continue
                
            logger.info(f"Processing {pdf_key}")
            
            # Extract text from PDF
            text = self.extract_text_from_pdf(pdf_key)
            
            # Chunk the text
            chunks = self.chunk_text(text)
            logger.info(f"Created {len(chunks)} chunks from {pdf_key}")
            
            # Map chunks to their source
            for chunk in chunks:
                chunk_to_source[chunk_idx] = {
                    "pdf": pdf_key,
                    "text": chunk
                }
                chunk_idx += 1
            
            all_chunks.extend(chunks)
            self.processed_files.add(pdf_key)
        
        # Generate embeddings - use distributed method for large chunk sets
        logger.info(f"Generating embeddings for {len(all_chunks)} chunks")
        if len(all_chunks) > 100:  # Threshold for distributed processing
            embeddings = self.get_embeddings_distributed(all_chunks)
        else:
            embeddings = self.get_embeddings(all_chunks)
        
        # Save the mapping between index and text
        self.index_to_text_map = chunk_to_source
        
        # Create FAISS index
        logger.info("Building FAISS index")
        embeddings_np = np.array(embeddings).astype('float32')
        dimension = embeddings_np.shape[1]
        index = IndexFlatL2(dimension)
        index.add(embeddings_np)
        
        # Store the index
        self.vector_store = {
            "index": index,
            "embeddings": embeddings_np
        }
        
        # Determine if we should use memory mapping for large indices
        if len(all_chunks) > LARGE_INDEX_THRESHOLD:
            self.use_mmap = True
            logger.info(f"Large index detected ({len(all_chunks)} vectors). Using memory mapping.")
        
        # Save to disk for persistence across Lambda invocations
        self.save_index_to_disk()
        
        # Save checkpoint to S3
        self.save_checkpoint()
        
        logger.info(f"Index built and saved with session ID: {self.session_id}")
    
    def process_large_pdfs(self, s3_keys: List[str]) -> None:
        """Process large PDF sets in batches with checkpointing."""
        total_pdfs = len(s3_keys)
        processed = 0
        
        # Check for existing checkpoint
        if self.load_checkpoint():
            logger.info(f"Loaded checkpoint with {len(self.processed_files)} already processed files")
            # Filter out already processed files
            s3_keys = [key for key in s3_keys if key not in self.processed_files]
            processed = total_pdfs - len(s3_keys)
            logger.info(f"Resuming from checkpoint: {processed}/{total_pdfs} files already processed")
        
        # Process in batches
        for i in range(0, len(s3_keys), BATCH_SIZE):
            batch = s3_keys[i:i+BATCH_SIZE]
            batch_size = len(batch)
            
            logger.info(f"Processing batch {i//BATCH_SIZE + 1}: {batch_size} files")
            
            all_chunks = []
            chunk_to_source = {}
            chunk_idx = len(self.index_to_text_map)  # Start from current index
            
            # 1. Extract text from each PDF in the batch
            for pdf_key in batch:
                logger.info(f"Processing {pdf_key}")
                
                # Extract text from PDF
                text = self.extract_text_from_pdf(pdf_key)
                
                # Chunk the text
                chunks = self.chunk_text(text)
                logger.info(f"Created {len(chunks)} chunks from {pdf_key}")
                
                # Map chunks to their source
                for chunk in chunks:
                    chunk_to_source[chunk_idx] = {
                        "pdf": pdf_key,
                        "text": chunk
                    }
                    chunk_idx += 1
                
                all_chunks.extend(chunks)
                self.processed_files.add(pdf_key)
            
            # 2. Generate embeddings for this batch
            if all_chunks:
                logger.info(f"Generating embeddings for {len(all_chunks)} chunks")
                if len(all_chunks) > 100:  # Threshold for distributed processing
                    embeddings = self.get_embeddings_distributed(all_chunks)
                else:
                    embeddings = self.get_embeddings(all_chunks)
                
                # 3. Update the index with new embeddings
                self.update_index(embeddings, chunk_to_source)
                
                # 4. Save checkpoint after each batch
                self.save_checkpoint()
            
            processed += batch_size
            logger.info(f"Processed {processed}/{total_pdfs} files")
            
            # Force garbage collection between batches
            gc.collect()
    
    def update_index(self, embeddings: List[List[float]], new_text_map: Dict[int, Dict]) -> None:
        """Update FAISS index incrementally with new embeddings."""
        # Convert embeddings to numpy array
        new_embeddings = np.array(embeddings).astype('float32')
        
        if self.vector_store is None:
            # Create new index if none exists
            dimension = new_embeddings.shape[1]
            index = IndexFlatL2(dimension)
            index.add(new_embeddings)
            
            self.vector_store = {
                "index": index,
                "embeddings": new_embeddings
            }
        else:
            # Merge with existing index
            current_embeddings = self.vector_store["embeddings"]
            combined = np.vstack([current_embeddings, new_embeddings])
            
            # Rebuild index with combined embeddings
            dimension = combined.shape[1]
            index = IndexFlatL2(dimension)
            index.add(combined)
            
            self.vector_store = {
                "index": index,
                "embeddings": combined
            }
        
        # Update text mapping
        self.index_to_text_map.update(new_text_map)
        
        # Determine if we should use memory mapping
        if self.vector_store["embeddings"].shape[0] > LARGE_INDEX_THRESHOLD:
            self.use_mmap = True
            logger.info(f"Large index detected ({self.vector_store['embeddings'].shape[0]} vectors). Using memory mapping.")
    
    def save_index_to_disk(self) -> None:
        """Save index and mapping to disk."""
        # Create directory if it doesn't exist
        os.makedirs(TEMP_DIR, exist_ok=True)
        
        # Save index
        if self.use_mmap:
            # For large indices, save the raw FAISS index
            index_file = os.path.join(TEMP_DIR, f"{self.session_id}_index.faiss")
            write_index(self.vector_store["index"], index_file)
            
            # Save embeddings separately using numpy's memory-mapped files
            embeddings_file = os.path.join(TEMP_DIR, f"{self.session_id}_embeddings.npy")
            np.save(embeddings_file, self.vector_store["embeddings"])
        else:
            # For smaller indices, use pickle
            vector_file = os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl")
            with open(vector_file, 'wb') as f:
                pickle.dump(self.vector_store, f)
        
        # Save mapping
        mapping_file = os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl")
        with open(mapping_file, 'wb') as f:
            pickle.dump(self.index_to_text_map, f)
    
    def save_checkpoint(self) -> None:
        """Save progress checkpoint to S3."""
        # Save metadata about the processing state
        checkpoint_data = {
            "session_id": self.session_id,
            "processed_files": list(self.processed_files),
            "last_activity": self.last_activity,
            "use_mmap": self.use_mmap,
            "vector_count": self.vector_store["embeddings"].shape[0] if self.vector_store else 0,
            "timestamp": datetime.now().isoformat()
        }
        
        # Save checkpoint metadata
        checkpoint_key = f"{CHECKPOINT_DIR}/{self.session_id}_metadata.json"
        s3_client.put_object(
            Bucket=self.bucket_name,
            Key=checkpoint_key,
            Body=json.dumps(checkpoint_data)
        )
        
        # Save the current index and mapping to S3
        if self.use_mmap:
            # Upload FAISS index
            index_file = os.path.join(TEMP_DIR, f"{self.session_id}_index.faiss")
            write_index(self.vector_store["index"], index_file)
            s3_client.upload_file(
                index_file, 
                self.bucket_name, 
                f"{CHECKPOINT_DIR}/{self.session_id}_index.faiss"
            )
            
            # Upload embeddings in chunks if they're large
            embeddings = self.vector_store["embeddings"]
            embeddings_file = os.path.join(TEMP_DIR, f"{self.session_id}_embeddings.npy")
            np.save(embeddings_file, embeddings)
            s3_client.upload_file(
                embeddings_file, 
                self.bucket_name, 
                f"{CHECKPOINT_DIR}/{self.session_id}_embeddings.npy"
            )
        else:
            # For smaller indices, use pickle
            vector_file = os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl")
            with open(vector_file, 'wb') as f:
                pickle.dump(self.vector_store, f)
            s3_client.upload_file(
                vector_file, 
                self.bucket_name, 
                f"{CHECKPOINT_DIR}/{self.session_id}_vectors.pkl"
            )
        
        # Upload mapping
        mapping_file = os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl")
        with open(mapping_file, 'wb') as f:
            pickle.dump(self.index_to_text_map, f)
        s3_client.upload_file(
            mapping_file, 
            self.bucket_name, 
            f"{CHECKPOINT_DIR}/{self.session_id}_mapping.pkl"
        )
        
        logger.info(f"Checkpoint saved to S3: {checkpoint_key}")
    
    def load_checkpoint(self) -> bool:
        """Load processing checkpoint from S3."""
        try:
            # Check if checkpoint exists
            checkpoint_key = f"{CHECKPOINT_DIR}/{self.session_id}_metadata.json"
            response = s3_client.get_object(
                Bucket=self.bucket_name,
                Key=checkpoint_key
            )
            
            # Load checkpoint metadata
            checkpoint_data = json.loads(response['Body'].read().decode('utf-8'))
            self.processed_files = set(checkpoint_data.get('processed_files', []))
            self.last_activity = checkpoint_data.get('last_activity', time.time())
            self.use_mmap = checkpoint_data.get('use_mmap', False)
            
            logger.info(f"Found checkpoint for session {self.session_id} with {len(self.processed_files)} processed files")
            
            # Load index and mapping
            if self.use_mmap:
                # Download FAISS index
                index_file = os.path.join(TEMP_DIR, f"{self.session_id}_index.faiss")
                s3_client.download_file(
                    self.bucket_name,
                    f"{CHECKPOINT_DIR}/{self.session_id}_index.faiss",
                    index_file
                )
                
                # Download embeddings
                embeddings_file = os.path.join(TEMP_DIR, f"{self.session_id}_embeddings.npy")
                s3_client.download_file(
                    self.bucket_name,
                    f"{CHECKPOINT_DIR}/{self.session_id}_embeddings.npy",
                    embeddings_file
                )
                
                # Load index and embeddings
                index = read_index(index_file)
                embeddings = np.load(embeddings_file)
                
                self.vector_store = {
                    "index": index,
                    "embeddings": embeddings
                }
            else:
                # Download vector store
                vector_file = os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl")
                s3_client.download_file(
                    self.bucket_name,
                    f"{CHECKPOINT_DIR}/{self.session_id}_vectors.pkl",
                    vector_file
                )
                
                # Load vector store
                with open(vector_file, 'rb') as f:
                    self.vector_store = pickle.load(f)
            
            # Download mapping
            mapping_file = os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl")
            s3_client.download_file(
                self.bucket_name,
                f"{CHECKPOINT_DIR}/{self.session_id}_mapping.pkl",
                mapping_file
            )
            
            # Load mapping
            with open(mapping_file, 'rb') as f:
                self.index_to_text_map = pickle.load(f)
            
            logger.info(f"Checkpoint loaded successfully with {len(self.index_to_text_map)} text chunks")
            return True
            
        except Exception as e:
            logger.info(f"No checkpoint found or error loading checkpoint: {str(e)}")
            return False
        
    def load_index(self, session_id: str) -> bool:
        """Load a previously created index and mapping."""
        # First check if there's a checkpoint
        self.session_id = session_id
        if self.load_checkpoint():
            return True
        
        # If no checkpoint, try the regular files
        vector_file = os.path.join(TEMP_DIR, f"{session_id}_vectors.pkl")
        mapping_file = os.path.join(TEMP_DIR, f"{session_id}_mapping.pkl")
        index_file = os.path.join(TEMP_DIR, f"{session_id}_index.faiss")
        embeddings_file = os.path.join(TEMP_DIR, f"{session_id}_embeddings.npy")
        
        # Check if we have memory-mapped files
        if os.path.exists(index_file) and os.path.exists(embeddings_file):
            logger.info(f"Loading memory-mapped index from local disk for session {session_id}")
            try:
                index = read_index(index_file)
                embeddings = np.load(embeddings_file)
                
                self.vector_store = {
                    "index": index,
                    "embeddings": embeddings
                }
                
                with open(mapping_file, 'rb') as f:
                    self.index_to_text_map = pickle.load(f)
                
                self.use_mmap = True
                return True
            except Exception as e:
                logger.error(f"Error loading memory-mapped index: {str(e)}")
        
        # Try to load regular pickle files from disk
        if os.path.exists(vector_file) and os.path.exists(mapping_file):
            logger.info(f"Loading index from local disk for session {session_id}")
            try:
                with open(vector_file, 'rb') as f:
                    self.vector_store = pickle.load(f)
                with open(mapping_file, 'rb') as f:
                    self.index_to_text_map = pickle.load(f)
                return True
            except Exception as e:
                logger.error(f"Error loading index from disk: {str(e)}")
        
        # If not on disk, try to load from S3
        try:
            logger.info(f"Loading index from S3 for session {session_id}")
            
            # First check for memory-mapped files
            try:
                s3_client.head_object(
                    Bucket=self.bucket_name,
                    Key=f"{CHECKPOINT_DIR}/{session_id}_index.faiss"
                )
                
                # Download FAISS index and embeddings
                s3_client.download_file(
                    self.bucket_name, 
                    f"{CHECKPOINT_DIR}/{session_id}_index.faiss", 
                    index_file
                )
                s3_client.download_file(
                    self.bucket_name, 
                    f"{CHECKPOINT_DIR}/{session_id}_embeddings.npy", 
                    embeddings_file
                )
                
                # Load index and embeddings
                index = read_index(index_file)
                embeddings = np.load(embeddings_file)
                
                self.vector_store = {
                    "index": index,
                    "embeddings": embeddings
                }
                
                # Download and load mapping
                s3_client.download_file(
                    self.bucket_name, 
                    f"{CHECKPOINT_DIR}/{session_id}_mapping.pkl", 
                    mapping_file
                )
                with open(mapping_file, 'rb') as f:
                    self.index_to_text_map = pickle.load(f)
                
                self.use_mmap = True
                return True
            except Exception:
                # Fall back to regular pickle files
                s3_vector_key = f"temp/{session_id}_vectors.pkl"
                s3_mapping_key = f"temp/{session_id}_mapping.pkl"
                
                s3_client.download_file(self.bucket_name, s3_vector_key, vector_file)
                s3_client.download_file(self.bucket_name, s3_mapping_key, mapping_file)
                
                with open(vector_file, 'rb') as f:
                    self.vector_store = pickle.load(f)
                with open(mapping_file, 'rb') as f:
                    self.index_to_text_map = pickle.load(f)
                
                return True
        except Exception as e:
            logger.error(f"Error loading index from S3: {str(e)}")
            return False
    
    def query(self, question: str, k: int = 5) -> Dict[str, Any]:
        """Query the RAG system with a natural language question.
        
        Args:
            question: User's question
            k: Number of relevant chunks to retrieve
            
        Returns:
            Dict with query results and answer
        """
        self.last_activity = time.time()
        
        # Free memory before query
        gc.collect()
        
        # Ensure vector store is loaded
        if not self.vector_store:
            raise ValueError("No index loaded. Build or load an index first.")
        
        # Generate embedding for the question
        question_embedding = self.get_embeddings([question])[0]
        question_embedding_np = np.array([question_embedding]).astype('float32')
        
        # Search for similar chunks
        try:
            D, I = self.vector_store["index"].search(question_embedding_np, k)
            
            # Get the relevant text chunks
            relevant_chunks = []
            for idx in I[0]:
                if idx in self.index_to_text_map:
                    relevant_chunks.append(self.index_to_text_map[idx])
            
            # Prepare context for the LLM
            context = "\n\n".join([chunk["text"] for chunk in relevant_chunks])
            
            # Create prompt for Claude
            prompt = f"""
            You are a helpful assistant that answers questions based on the provided context.
            
            Context:
            {context}
            
            Question: {question}
            
            Answer the question based on the context provided. If the answer is not in the context, say "I don't have enough information to answer this question."
            Include citations to the source documents where appropriate using [Source: document_name].
            """
            
            # Call Claude via AWS Bedrock
            try:
                request_body = json.dumps({
                    "anthropic_version": "bedrock-2023-05-31",
                    "max_tokens": 1000,
                    "messages": [
                        {"role": "user", "content": prompt}
                    ]
                })
                
                response = bedrock_client.invoke_model(
                    modelId=CLAUDE_MODEL_ID,
                    contentType="application/json",
                    accept="application/json",
                    body=request_body
                )
                
                response_body = json.loads(response['body'].read())
                answer = response_body['content'][0]['text']
                
                # Update checkpoint after successful query
                if hasattr(self, 'processed_files') and self.processed_files:
                    self.save_checkpoint()
                
                return {
                    "question": question,
                    "answer": answer,
                    "sources": [chunk["pdf"] for chunk in relevant_chunks],
                    "context": context
                }
                
            except Exception as e:
                logger.error(f"Error calling Bedrock: {str(e)}")
                return {
                    "question": question,
                    "answer": "Sorry, I encountered an error while processing your question.",
                    "error": str(e)
                }
        except Exception as e:
            logger.error(f"Error during vector search: {str(e)}")
            return {
                "question": question,
                "answer": "Sorry, I encountered an error during the search process.",
                "error": str(e)
            }
    
    def clean_up(self) -> None:
        """Clean up temporary resources."""
        logger.info(f"Cleaning up session {self.session_id}")
        
        # Remove local files
        files_to_remove = [
            os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl"),
            os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl"),
            os.path.join(TEMP_DIR, f"{self.session_id}_index.faiss"),
            os.path.join(TEMP_DIR, f"{self.session_id}_embeddings.npy")
        ]
        
        for file_path in files_to_remove:
            if os.path.exists(file_path):
                try:
                    os.remove(file_path)
                except Exception as e:
                    logger.warning(f"Error removing {file_path}: {str(e)}")
        
        # Remove S3 temporary files and checkpoints
        s3_keys_to_delete = [
            f"temp/{self.session_id}_vectors.pkl",
            f"temp/{self.session_id}_mapping.pkl",
            f"{CHECKPOINT_DIR}/{self.session_id}_metadata.json",
            f"{CHECKPOINT_DIR}/{self.session_id}_index.faiss",
            f"{CHECKPOINT_DIR}/{self.session_id}_embeddings.npy",
            f"{CHECKPOINT_DIR}/{self.session_id}_mapping.pkl",
            f"{CHECKPOINT_DIR}/{self.session_id}_vectors.pkl"
        ]
        
        for key in s3_keys_to_delete:
            try:
                s3_client.delete_object(
                    Bucket=self.bucket_name,
                    Key=key
                )
            except Exception as e:
                logger.warning(f"Error deleting S3 object {key}: {str(e)}")
        
        logger.info("Cleanup completed")
    
    def check_session_timeout(self) -> bool:
        """Check if the session has timed out."""
        current_time = time.time()
        elapsed = current_time - self.last_activity
        return elapsed > SESSION_TIMEOUT


# Lambda handler for embedding generation
def embedding_handler(event, context):
    """Lambda handler for distributed embedding generation."""
    try:
        texts = event.get('texts', [])
        model_id = event.get('model_id', EMBEDDING_MODEL_ID)
        
        if not texts:
            return {
                'statusCode': 400,
                'error': 'No texts provided'
            }
        
        logger.info(f"Generating embeddings for {len(texts)} chunks")
        
        embeddings = []
        for text in texts:
            try:
                # Format the request for the embedding model
                request_body = json.dumps({
                    "inputText": text
                })
                
                response = bedrock_client.invoke_model(
                    modelId=model_id,
                    contentType="application/json",
                    accept="application/json",
                    body=request_body
                )
                
                response_body = json.loads(response['body'].read())
                embedding_vector = response_body['embedding']
                embeddings.append(embedding_vector)
                
            except Exception as e:
                logger.error(f"Error generating embedding: {str(e)}")
                # Return a zero vector as fallback
                embeddings.append([0.0] * 1536)  # Default size for many embedding models
        
        return {
            'statusCode': 200,
            'embeddings': embeddings
        }
    except Exception as e:
        logger.error(f"Error in embedding handler: {str(e)}")
        return {
            'statusCode': 500,
            'error': str(e)
        }

# Lambda handler functions for serverless operation
def upload_handler(event, context):
    """Lambda handler for PDF upload."""
    try:
        bucket = event['Records'][0]['s3']['bucket']['name']
        key = event['Records'][0]['s3']['object']['key']
        
        # Only process PDFs
        if not key.lower().endswith('.pdf'):
            return {
                'statusCode': 200,
                'body': json.dumps('Not a PDF, skipping')
            }
        
        # Create a new RAG instance
        rag = EnhancedServerlessRAG(bucket_name=bucket)
        
        # Process the uploaded PDF
        rag.build_index([key])
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': f'Successfully processed {key}',
                'session_id': rag.session_id
            })
        }
    except Exception as e:
        logger.error(f"Error in upload handler: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e)
            })
        }

def query_handler(event, context):
    """Lambda handler for querying the RAG system."""
    try:
        body = json.loads(event['body'])
        session_id = body.get('session_id')
        question = body.get('question')
        bucket_name = body.get('bucket', DEFAULT_BUCKET_NAME)
        
        if not session_id or not question:
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'error': 'Missing required parameters: session_id and question'
                })
            }
        
        # Create a RAG instance and load the specified session
        rag = EnhancedServerlessRAG(bucket_name=bucket_name, session_id=session_id)
        if not rag.load_index(session_id):
            return {
                'statusCode': 404,
                'body': json.dumps({
                    'error': f'Session {session_id} not found'
                })
            }
        
        # Process the query
        result = rag.query(question)
        
        return {
            'statusCode': 200,
            'body': json.dumps(result)
        }
    except Exception as e:
        logger.error(f"Error in query handler: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e)
            })
        }

def cleanup_handler(event, context):
    """Lambda handler for cleaning up resources."""
    try:
        body = json.loads(event['body'])
        session_id = body.get('session_id')
        bucket_name = body.get('bucket', DEFAULT_BUCKET_NAME)
        
        if not session_id:
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'error': 'Missing required parameter: session_id'
                })
            }
        
        # Create a RAG instance and load the specified session
        rag = EnhancedServerlessRAG(bucket_name=bucket_name, session_id=session_id)
        if rag.load_index(session_id):
            rag.clean_up()
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': f'Successfully cleaned up session {session_id}'
            })
        }
    except Exception as e:
        logger.error(f"Error in cleanup handler: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e)
            })
        }

def scheduled_cleanup_handler(event, context):
    """Lambda handler for scheduled cleanup of expired sessions."""
    try:
        bucket_name = event.get('bucket', DEFAULT_BUCKET_NAME)
        
        # List all temp files and checkpoints in S3
        prefixes_to_check = ['temp/', f'{CHECKPOINT_DIR}/']
        deleted_files = 0
        
        for prefix in prefixes_to_check:
            response = s3_client.list_objects_v2(
                Bucket=bucket_name,
                Prefix=prefix
            )
            
            if 'Contents' not in response:
                continue
            
            # Get the current time
            current_time = datetime.now()
            
            # Check each file's last modified time
            for obj in response['Contents']:
                key = obj['Key']
                last_modified = obj['LastModified']
                
                # If the file is older than the session timeout, delete it
                if current_time - last_modified > timedelta(seconds=SESSION_TIMEOUT):
                    s3_client.delete_object(
                        Bucket=bucket_name,
                        Key=key
                    )
                    logger.info(f"Deleted expired file: {key}")
                    deleted_files += 1
        
        # Also check for checkpoint metadata files to find sessions with auto-delete
        try:
            response = s3_client.list_objects_v2(
                Bucket=bucket_name,
                Prefix=f'{CHECKPOINT_DIR}/'
            )
            
            if 'Contents' in response:
                for obj in response['Contents']:
                    key = obj['Key']
                    
                    # Only process metadata files
                    if not key.endswith('_metadata.json'):
                        continue
                    
                    # Get the metadata
                    try:
                        metadata_obj = s3_client.get_object(
                            Bucket=bucket_name,
                            Key=key
                        )
                        metadata = json.loads(metadata_obj['Body'].read().decode('utf-8'))
                        
                        # Check if this session has been inactive
                        last_activity = metadata.get('last_activity', 0)
                        if time.time() - last_activity > SESSION_TIMEOUT:
                            # Extract session ID from the key
                            session_id = key.split('/')[-1].replace('_metadata.json', '')
                            
                            # Create a RAG instance to clean up this session
                            rag = EnhancedServerlessRAG(bucket_name=bucket_name, session_id=session_id)
                            rag.clean_up()
                            logger.info(f"Cleaned up inactive session: {session_id}")
                    except Exception as e:
                        logger.warning(f"Error processing metadata file {key}: {str(e)}")
        except Exception as e:
            logger.warning(f"Error checking checkpoint metadata: {str(e)}")
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': f'Scheduled cleanup completed. Deleted {deleted_files} expired files.'
            })
        }
    except Exception as e:
        logger.error(f"Error in scheduled cleanup handler: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e)
            })
        }
# Batch processing handler for large PDF sets
def batch_processing_handler(event, context):
    """Lambda handler for processing large batches of PDFs."""
    try:
        bucket_name = event.get('bucket', DEFAULT_BUCKET_NAME)
        session_id = event.get('session_id')
        pdf_keys = event.get('pdf_keys', [])
        
        if not pdf_keys:
            return {
                'statusCode': 400,
                'body': json.dumps({
                    'error': 'No PDF keys provided'
                })
            }
        
        # Create a RAG instance
        rag = EnhancedServerlessRAG(bucket_name=bucket_name, session_id=session_id)
        
        # Process the PDFs in batches
        rag.process_large_pdfs(pdf_keys)
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': f'Successfully processed {len(pdf_keys)} PDFs',
                'session_id': rag.session_id
            })
        }
    except Exception as e:
        logger.error(f"Error in batch processing handler: {str(e)}")
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': str(e)
            })
        }

# Example usage script (not part of Lambda)
if __name__ == "__main__":
    # Initialize
    rag = EnhancedServerlessRAG()
    
    # Get file path from command line
    import sys
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        print("Usage: python3 bedrock-lite-rag.py <path-to-file-or-directory>")
        sys.exit(1)
    
    # Upload and process
    uploaded_keys = rag.upload_pdf(path)
    if uploaded_keys:
        rag.build_index(uploaded_keys)
        
        # Interactive query
        while True:
            question = input("\nEnter question (or 'exit' to quit): ")
            if question.lower() == 'exit':
                break
            answer = rag.query(question)
            print(f"\nAnswer: {answer['answer']}")
            if answer['sources']:
                print("\nSources:")
                for source in list(set(answer['sources'])):  # Deduplicate
                    print(f"- {source}")
        
        rag.clean_up()
    else:
        print("No PDF files found to process")
