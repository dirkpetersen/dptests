from flask import Flask, render_template, request, make_response
import uuid
import boto3
from markitdown import MarkItDown
import io
import os
import json
import logging
import hashlib
from typing import Tuple, List, Optional
from pathlib import Path
from werkzeug.utils import secure_filename
from botocore.exceptions import BotoCoreError, ClientError
from botocore.config import Config
from datetime import datetime, timezone, timedelta

from config import *

# Configure logging based on Flask debug mode
logger = logging.getLogger(__name__)
def configure_logging():
    log_level = logging.DEBUG if app.debug else logging.INFO
    logging.basicConfig(level=log_level, format=LOG_FORMAT)

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

def get_user_id():
    """Get or create user ID from cookie"""
    user_id = request.cookies.get(COOKIE_NAME)
    if not user_id:
        user_id = str(uuid.uuid4())
    return user_id

def get_policy_cache_path(user_id: str) -> Path:
    """Get path to user's cached policy file"""
    return Path(CACHE_DIR) / f"policy_{user_id}.txt"

def save_policy_to_cache(text: str, user_id: str):
    """Save policy text to user-specific cache"""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache_path = get_policy_cache_path(user_id)
    cache_path.write_text(text)
    logger.debug(f"Saved policy to cache: {cache_path}")

def get_cached_policy(user_id: str) -> Optional[str]:
    """Get cached policy text for user"""
    cache_path = get_policy_cache_path(user_id)
    if cache_path.exists():
        logger.debug(f"Policy cache hit: {cache_path}")
        return cache_path.read_text()
    return None

def extract_text_from_pdf(pdf_file) -> str:
    """
    Extract text content from a PDF file using MarkItDown.
    
    Args:
        pdf_file: File object containing PDF data
        
    Returns:
        str: Extracted text from the PDF in markdown format
        
    Raises:
        Exception: If PDF parsing fails
    """
    # Save the uploaded file temporarily
    temp_path = "temp.pdf"
    pdf_file.save(temp_path)
    
    try:
        md = MarkItDown()
        result = md.convert(temp_path)
        return result.text_content
    finally:
        # Clean up temporary file
        if os.path.exists(temp_path):
            os.remove(temp_path)

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
    # Create simple analysis prompt with full documents
    analysis_prompt = f"""Human: Compare these two documents and determine if the second document meets the requirements specified in the first document.

Policy Document:
{policy_text}

Submission Document:
{submission_text}

Based on all these comparisons, respond with exactly one word (GREEN, YELLOW, ORANGE, or RED).
In addition provide an explanation on how specific requirements are met (GREEN) or may not be met.

GREEN means all requirements are fully met.
YELLOW means all quantifiable/numerical requirements are met but other requirements are ambiguous.
ORANGE means both numerical and other requirements are ambiguous and need clarification.
RED means one or more requirements are clearly not met.

""" 
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
        
        if status in ["GREEN", "YELLOW", "ORANGE", "RED"]:
            return status, explanation
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
            user_id = get_user_id()
            
            # Handle policy document
            if policy_file and policy_file.filename:
                policy_text = extract_text_from_pdf(policy_file)
                save_policy_to_cache(policy_text, user_id)
            else:
                # Try to load cached policy
                policy_text = get_cached_policy(user_id)
                if not policy_text:
                    return render_template('index.html', error="No policy document available. Please upload one.")
                logger.debug("Using cached policy document")
            
            submission_text = extract_text_from_pdf(submission_file)
            
            result, explanation = evaluate_requirements(policy_text, submission_text)
            logger.info(f"Evaluation result: {result}, Explanation: {explanation}")
            
            # Create response with cookie
            response = make_response(render_template('index.html', 
                                                   result=result, 
                                                   explanation=explanation))
            response.set_cookie(COOKIE_NAME, user_id, max_age=COOKIE_MAX_AGE)
            return response
            
        except Exception as e:
            logger.error(f"Error processing request: {str(e)}")
            return render_template('index.html', error=str(e))
    
    return render_template('index.html')

if __name__ == '__main__':
    configure_logging()
    app.run(host='0.0.0.0', port=5000, debug=True)
