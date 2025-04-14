"""
Serverless RAG Application using AWS Bedrock and S3
===================================================

This implementation provides a complete solution for:
1. Setting up the required AWS resources
2. Processing PDFs from S3
3. Creating temporary vector stores
4. Querying the knowledge base with AWS Bedrock
5. Tearing down resources when done

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
from typing import List, Dict, Any, Optional
import tempfile
import shutil
from datetime import datetime, timedelta

# For PDF processing
import PyPDF2
from langchain.text_splitter import RecursiveCharacterTextSplitter

# For vector operations
import numpy as np
from faiss import IndexFlatL2
import pickle

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
SESSION_TIMEOUT = 60 * 60  # 1 hour in seconds
CLAUDE_MODEL_ID = "anthropic.claude-3-sonnet-20240229-v1:0"  # Update as needed
EMBEDDING_MODEL_ID = "amazon.titan-embed-text-v1"  # Update as needed
CHUNK_SIZE = 1000
CHUNK_OVERLAP = 200

class ServerlessRAG:
    def __init__(self, bucket_name: str = DEFAULT_BUCKET_NAME):
        """Initialize the RAG application."""
        self.bucket_name = bucket_name
        self.session_id = str(uuid.uuid4())
        self.vector_store = None
        self.index_to_text_map = {}
        self.last_activity = time.time()
        
        # Create temporary directory for storing vector data
        os.makedirs(TEMP_DIR, exist_ok=True)
        
        # Check if bucket exists, if not create it
        self._ensure_bucket_exists()
    
    def _ensure_bucket_exists(self):
        """Create S3 bucket if it doesn't exist."""
        try:
            s3_client.head_bucket(Bucket=self.bucket_name)
            logger.info(f"Bucket {self.bucket_name} already exists")
        except:
            logger.info(f"Creating bucket {self.bucket_name}")
            s3_client.create_bucket(Bucket=self.bucket_name)
            
            # Enable versioning for data protection
            s3_client.put_bucket_versioning(
                Bucket=self.bucket_name,
                VersioningConfiguration={'Status': 'Enabled'}
            )
    
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
                
                for page_num in range(len(pdf_reader.pages)):
                    page = pdf_reader.pages[page_num]
                    text += page.extract_text() + "\n\n"
            
            return text
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
        
        for i, text in enumerate(texts):
            if i % 10 == 0:
                logger.info(f"Generating embedding for chunk {i+1}/{len(texts)}")
            
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
                embeddings.append(embedding_vector)
                
            except Exception as e:
                logger.error(f"Error generating embedding: {str(e)}")
                # Return a zero vector as fallback
                embeddings.append([0.0] * 1536)  # Default size for many embedding models
        
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
        
        all_chunks = []
        chunk_to_source = {}
        chunk_idx = 0
        
        # Process each PDF
        for pdf_key in pdfs:
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
        
        # Generate embeddings
        logger.info(f"Generating embeddings for {len(all_chunks)} chunks")
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
        
        # Save to disk for persistence across Lambda invocations
        vector_file = os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl")
        with open(vector_file, 'wb') as f:
            pickle.dump(self.vector_store, f)
        
        mapping_file = os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl")
        with open(mapping_file, 'wb') as f:
            pickle.dump(self.index_to_text_map, f)
        
        # Optional: Upload to S3 for longer persistence
        s3_client.upload_file(vector_file, self.bucket_name, f"temp/{self.session_id}_vectors.pkl")
        s3_client.upload_file(mapping_file, self.bucket_name, f"temp/{self.session_id}_mapping.pkl")
        
        logger.info(f"Index built and saved with session ID: {self.session_id}")
        
    def load_index(self, session_id: str) -> bool:
        """Load a previously created index and mapping."""
        vector_file = os.path.join(TEMP_DIR, f"{session_id}_vectors.pkl")
        mapping_file = os.path.join(TEMP_DIR, f"{session_id}_mapping.pkl")
        
        # Try to load from disk first
        if os.path.exists(vector_file) and os.path.exists(mapping_file):
            logger.info(f"Loading index from local disk for session {session_id}")
            with open(vector_file, 'rb') as f:
                self.vector_store = pickle.load(f)
            with open(mapping_file, 'rb') as f:
                self.index_to_text_map = pickle.load(f)
            self.session_id = session_id
            return True
        
        # If not on disk, try to load from S3
        try:
            logger.info(f"Loading index from S3 for session {session_id}")
            s3_vector_key = f"temp/{session_id}_vectors.pkl"
            s3_mapping_key = f"temp/{session_id}_mapping.pkl"
            
            s3_client.download_file(self.bucket_name, s3_vector_key, vector_file)
            s3_client.download_file(self.bucket_name, s3_mapping_key, mapping_file)
            
            with open(vector_file, 'rb') as f:
                self.vector_store = pickle.load(f)
            with open(mapping_file, 'rb') as f:
                self.index_to_text_map = pickle.load(f)
            
            self.session_id = session_id
            return True
        except Exception as e:
            logger.error(f"Error loading index: {str(e)}")
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
        
        # Ensure vector store is loaded
        if not self.vector_store:
            raise ValueError("No index loaded. Build or load an index first.")
        
        # Generate embedding for the question
        question_embedding = self.get_embeddings([question])[0]
        question_embedding_np = np.array([question_embedding]).astype('float32')
        
        # Search for similar chunks
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
    
    def clean_up(self) -> None:
        """Clean up temporary resources."""
        logger.info(f"Cleaning up session {self.session_id}")
        
        # Remove local files
        vector_file = os.path.join(TEMP_DIR, f"{self.session_id}_vectors.pkl")
        mapping_file = os.path.join(TEMP_DIR, f"{self.session_id}_mapping.pkl")
        
        if os.path.exists(vector_file):
            os.remove(vector_file)
        if os.path.exists(mapping_file):
            os.remove(mapping_file)
        
        # Remove S3 temporary files
        try:
            s3_client.delete_object(
                Bucket=self.bucket_name,
                Key=f"temp/{self.session_id}_vectors.pkl"
            )
            s3_client.delete_object(
                Bucket=self.bucket_name,
                Key=f"temp/{self.session_id}_mapping.pkl"
            )
        except Exception as e:
            logger.warning(f"Error cleaning up S3: {str(e)}")
        
        logger.info("Cleanup completed")
    
    def check_session_timeout(self) -> bool:
        """Check if the session has timed out."""
        current_time = time.time()
        elapsed = current_time - self.last_activity
        return elapsed > SESSION_TIMEOUT


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
        rag = ServerlessRAG(bucket_name=bucket)
        
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
        rag = ServerlessRAG(bucket_name=bucket_name)
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
        rag = ServerlessRAG(bucket_name=bucket_name)
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
        
        # List all temp files in S3
        response = s3_client.list_objects_v2(
            Bucket=bucket_name,
            Prefix='temp/'
        )
        
        if 'Contents' not in response:
            return {
                'statusCode': 200,
                'body': json.dumps({
                    'message': 'No temporary files to clean up'
                })
            }
        
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
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': 'Scheduled cleanup completed'
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


# Example usage script (not part of Lambda)
if __name__ == "__main__":
    # Example: Manual local usage
    
    # Initialize
    rag = ServerlessRAG()
    
    # Upload a PDF (if running locally)
    pdf_key = rag.upload_pdf("sample.pdf")
    
    # Process PDFs and build index
    rag.build_index([pdf_key])
    
    # Query the system
    answer = rag.query("What is the main topic of this document?")
    print(json.dumps(answer, indent=2))
    
    # Clean up when done
    rag.clean_up()
