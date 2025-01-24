from flask import Flask, render_template, request
import boto3
import PyPDF2
import io
import os
import json
import logging
import hashlib
from typing import Tuple, List, Optional
from pathlib import Path
from werkzeug.utils import secure_filename
from langchain_community.embeddings import BedrockEmbeddings
from langchain.text_splitter import RecursiveCharacterTextSplitter
from langchain_community.vectorstores import FAISS
from botocore.exceptions import BotoCoreError, ClientError
from botocore.config import Config

from config import *

# Configure logging
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH
app.config['SECRET_KEY'] = os.urandom(24)

def get_bedrock_client():
    """Initialize Bedrock client with local credentials and retry configuration"""
    try:
        session = boto3.Session()
        region = session.region_name or BEDROCK_REGION
        
        # Configure retry behavior
        retry_config = Config(
            retries={
                "max_attempts": MAX_RETRIES,
                "mode": "standard",
            }
        )
        
        # Create Bedrock client using local credentials
        client = session.client(
            service_name='bedrock-runtime',
            region_name=region,
            config=retry_config
        )
        
        return client
    except (BotoCoreError, ClientError) as e:
        app.logger.error(f"Failed to initialize Bedrock client: {str(e)}")
        raise

# Initialize Bedrock client
bedrock = get_bedrock_client()

# Configure embeddings and text splitter
embeddings = BedrockEmbeddings(
    client=bedrock,
    model_id="amazon.titan-embed-text-v1"
)

text_splitter = RecursiveCharacterTextSplitter(
    chunk_size=CHUNK_SIZE,
    chunk_overlap=CHUNK_OVERLAP
)

def get_cache_path(file_content: bytes) -> Path:
    """Generate cache path for file content"""
    file_hash = hashlib.md5(file_content).hexdigest()
    return Path(CACHE_DIR) / f"{file_hash}.txt"

def get_cached_text(file_content: bytes) -> Optional[str]:
    """Get cached text content if available"""
    if not CACHE_ENABLED:
        return None
        
    cache_path = get_cache_path(file_content)
    if cache_path.exists():
        logger.debug(f"Cache hit: {cache_path}")
        return cache_path.read_text()
    return None

def save_to_cache(file_content: bytes, text: str, is_policy: bool = False):
    """Save extracted text to cache"""
    if not CACHE_ENABLED:
        return
        
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = get_cache_path(file_content)
    cache_path.write_text(text)
    logger.debug(f"Saved to cache: {cache_path}")
    
    # If it's a policy document, also save as the latest policy
    if is_policy:
        policy_path = Path(CACHE_DIR) / POLICY_CACHE_FILE
        policy_path.write_text(text)
        logger.debug(f"Saved latest policy: {policy_path}")

def extract_text_from_pdf(pdf_file) -> str:
    """
    Extract text content from a PDF file.
    
    Args:
        pdf_file: File object containing PDF data
        
    Returns:
        str: Extracted text from the PDF
        
    Raises:
        PyPDF2.PdfReadError: If PDF parsing fails
    """
    # Read file content for caching
    file_content = pdf_file.read()
    pdf_file.seek(0)  # Reset file pointer
    
    # Check cache
    cached_text = get_cached_text(file_content)
    if cached_text is not None:
        return cached_text
    
    # Extract text if not cached
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text()
    
    # Save to cache
    save_to_cache(file_content, text)
    return text

def evaluate_requirements(policy_text: str, submission_text: str) -> Tuple[str, str]:
    """
    Evaluate if a submission document meets the requirements in a policy document.
    
    Args:
        policy_text: Text content of the policy document
        submission_text: Text content of the submission document
        
    Returns:
        Tuple[str, str]: (result status, explanation if status is YELLOW)
        
    Raises:
        BotoCoreError: If AWS Bedrock API call fails
        ValueError: If response parsing fails
    """
    # Split documents into chunks
    policy_chunks = text_splitter.split_text(policy_text)
    submission_chunks = text_splitter.split_text(submission_text)
    
    # Create vector store from policy document
    policy_store = FAISS.from_texts(policy_chunks, embeddings)
    
    # For each submission chunk, find relevant policy requirements
    relevant_pairs = []
    for chunk in submission_chunks:
        similar_docs = policy_store.similarity_search(chunk, k=2)
        relevant_pairs.extend([(chunk, doc.page_content) for doc in similar_docs])
    
    # Create analysis prompt with relevant chunks
    analysis_prompt = """Human: I will provide pairs of text chunks from two documents. 
For each pair, the first is from a submission document and the second is from a policy document. 
Analyze if the submission meets the requirements in the policy.\n\n"""
    
    for sub_chunk, pol_chunk in relevant_pairs:
        analysis_prompt += f"Submission chunk:\n{sub_chunk}\n\nMatching policy chunk:\n{pol_chunk}\n\n"
    
    analysis_prompt += "Based on all these comparisons, respond with exactly one word (GREEN, YELLOW, or RED) "
    analysis_prompt += "followed by a brief explanation if YELLOW. GREEN means the submission fully meets requirements, "
    analysis_prompt += "RED means it clearly doesn't, and YELLOW means there are uncertainties that need human review."
    
    try:
        request_body = {
            "anthropic_version": ANTHROPIC_VERSION,
            "max_tokens": MAX_TOKENS,
            "temperature": TEMPERATURE,
            "messages": [
                {
                    "role": "user",
                    "content": analysis_prompt
                }
            ]
        }
        logger.debug(f"Bedrock request body: {json.dumps(request_body, indent=2)}")
        
        response = bedrock.invoke_model(
            modelId=MODEL_ID,
            body=json.dumps(request_body)
        )
        
        response_body = json.loads(response['body'].read())
        logger.debug(f"Bedrock response: {json.dumps(response_body, indent=2)}")
        
        # Extract text content from Claude 3's response structure
        content = response_body.get('content', [])
        if not content or not isinstance(content, list):
            logger.error(f"Unexpected response structure: {response_body}")
            raise ValueError("Response missing content array")
            
        # Get the text from the first content item
        text_content = next((item['text'] for item in content if item['type'] == 'text'), None)
        if not text_content:
            logger.error(f"No text content found in response: {response_body}")
            raise ValueError("Response missing text content")
            
        # Process the response text
        # First word of first line is the status
        lines = text_content.strip().split('\n')
        first_line = lines[0].strip()
        status = first_line.split('.')[0].strip().upper()
        
        # Combine the rest of first line (after status) with remaining lines for explanation
        remaining_first_line = first_line.split('.', 1)[1].strip() if '.' in first_line else ""
        remaining_lines = lines[1:] if len(lines) > 1 else []
        
        explanation = remaining_first_line
        if remaining_lines:
            explanation += "\n" + '\n'.join(remaining_lines)
        
        if status in ["GREEN", "YELLOW", "RED"]:
            return status, explanation if status == "YELLOW" else ""
        else:
            logger.error(f"Unexpected status in response: {text_content}")
            raise ValueError(f"Invalid status in response: {status}")
            
    except (BotoCoreError, ClientError) as e:
        logger.error(f"Bedrock API error: {str(e)}")
        raise
    except (KeyError, json.JSONDecodeError, ValueError) as e:
        logger.error(f"Failed to parse Bedrock response: {str(e)}")
        raise

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if 'submission' not in request.files:
            return render_template('index.html', error="Submission document is required")
        
        submission_file = request.files['submission']
        if submission_file.filename == '':
            return render_template('index.html', error="Submission document is required")
            
        policy_file = request.files.get('policy')
        
        try:
            # Handle policy document
            if policy_file.filename:
                policy_text = extract_text_from_pdf(policy_file)
                save_to_cache(policy_file.read(), policy_text, is_policy=True)
            else:
                # Try to load the latest policy
                policy_path = Path(CACHE_DIR) / POLICY_CACHE_FILE
                if not policy_path.exists():
                    return render_template('index.html', error="No policy document available. Please upload one.")
                policy_text = policy_path.read_text()
                logger.debug("Using cached policy document")
            
            submission_text = extract_text_from_pdf(submission_file)
            
            result, explanation = evaluate_requirements(policy_text, submission_text)
            logger.info(f"Evaluation result: {result}, Explanation: {explanation}")
            return render_template('index.html', 
                                result=result, 
                                explanation=explanation if result == "YELLOW" else None)
            
        except Exception as e:
            logger.error(f"Error processing request: {str(e)}")
            return render_template('index.html', error=str(e))
    
    return render_template('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
