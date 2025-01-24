from flask import Flask, render_template, request
import boto3
import PyPDF2
import io
import os
import json
import logging
from typing import Tuple, List
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
    pdf_reader = PyPDF2.PdfReader(pdf_file)
    text = ""
    for page in pdf_reader.pages:
        text += page.extract_text()
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
        text_upper = text_content.upper()
        if "GREEN" in text_upper:
            return "GREEN", ""
        elif "RED" in text_upper:
            return "RED", ""
        else:
            # Extract explanation after "YELLOW"
            parts = text_content.split("YELLOW", 1)
            explanation = parts[1].strip() if len(parts) > 1 else text_content
        return "YELLOW", explanation

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if 'policy' not in request.files or 'submission' not in request.files:
            return render_template('index.html', error="Both files are required")
        
        policy_file = request.files['policy']
        submission_file = request.files['submission']
        
        if policy_file.filename == '' or submission_file.filename == '':
            return render_template('index.html', error="Both files are required")
        
        try:
            policy_text = extract_text_from_pdf(policy_file)
            submission_text = extract_text_from_pdf(submission_file)
            
            result, explanation = evaluate_requirements(policy_text, submission_text)
            return render_template('index.html', result=result, explanation=explanation)
            
        except Exception as e:
            return render_template('index.html', error=str(e))
    
    return render_template('index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
