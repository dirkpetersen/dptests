#! /usr/bin/env python3

import os, json, uuid, tempfile, shutil, atexit, argparse, random, string,  logging, textwrap
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from flask import Flask, render_template, request, jsonify, session, redirect, url_for
from flask_session import Session
from werkzeug.utils import secure_filename
import boto3
from botocore.exceptions import ClientError
from botocore.config import Config
import markdown
import fitz  # PyMuPDF
import re
from dotenv import load_dotenv
import signal
import sys
import threading
from threading import Lock

# Load environment variables from .env file
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Parse command line arguments
parser = argparse.ArgumentParser(description='BedBot - AI Chat Assistant')
parser.add_argument('--no-bucket', action='store_true', help='Use local filesystem instead of S3 bucket')
parser.add_argument('--debug', action='store_true', help='Enable debug mode to print API messages')
parser.add_argument('--model', action='store', help='Bedrock model to use', default='us.amazon.nova-premier-v1:0')
args = parser.parse_args()

# Configuration
USE_S3_BUCKET = not args.no_bucket
DEBUG_MODE = args.debug
MAX_FILE_SIZE = 4.5 * 1024 * 1024  # 4.5 MB per file
MAX_FILES_PER_SESSION = 1000
BEDROCK_MODEL = os.getenv('BEDROCK_MODEL', args.model)
# PDF merging configuration from environment
PDF_MERGE_ENABLED = os.getenv('PDF_MERGE', '0').strip().lower() in ('1', 'true', 'yes', 'on')
# PDF conversion configuration - temporarily disable to test
PDF_CONVERSION_ENABLED = os.getenv('PDF_CONVERSION', '1').strip().lower() in ('1', 'true', 'yes', 'on')
# Maximum concurrent Bedrock connections for parallel PDF conversion
BEDROCK_MAX_CONCURRENT = int(os.getenv('BEDROCK_MAXCON', '3'))
# PDF conversion model - use different model for PDF conversion if specified
PDF_CONVERT_MODEL = os.getenv('PDF_CONVERT_MODEL', BEDROCK_MODEL)

if DEBUG_MODE:
    logger.info(f"PDF merging enabled: {PDF_MERGE_ENABLED} (PDF_MERGE={os.getenv('PDF_MERGE', '0')})")
    logger.info(f"PDF conversion enabled: {PDF_CONVERSION_ENABLED} (PDF_CONVERSION={os.getenv('PDF_CONVERSION', '1')})")
    logger.info(f"Max concurrent Bedrock connections: {BEDROCK_MAX_CONCURRENT} (BEDROCK_MAXCON={os.getenv('BEDROCK_MAXCON', '3')})")
    logger.info(f"PDF conversion model: {PDF_CONVERT_MODEL} (PDF_CONVERT_MODEL={os.getenv('PDF_CONVERT_MODEL', 'default')})")
    logger.info(f"Chat model: {BEDROCK_MODEL} (BEDROCK_MODEL={os.getenv('BEDROCK_MODEL', 'default')}, --model={args.model})")
log_args_info = f"Using command line args: --no-bucket={args.no_bucket}, --debug={args.debug}, --model={args.model}"

app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'your-secret-key-change-this')
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024  # 1GB max file size

# Configure Flask-Session for server-side sessions
app.config['SESSION_TYPE'] = 'filesystem'
app.config['SESSION_PERMANENT'] = False
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_KEY_PREFIX'] = 'bedbot:'
app.config['SESSION_FILE_THRESHOLD'] = 100  # Max number of sessions before cleanup

# Flask session directory will be created on demand
temp_session_dir = None

# Track all session upload folders for cleanup (local mode only)
session_upload_folders = set()

# Remove complex locking - use Flask's built-in session handling

# S3 bucket for file storage (S3 mode only)
s3_bucket_name = None
bucket_name_file = None

# Cleanup flag to prevent duplicate cleanup
cleanup_performed = False

# Initialize Flask-Session (will create session dir when first needed)
Session(app)

def generate_bucket_name():
    """Generate a unique bucket name with random suffix"""
    random_suffix = ''.join(random.choices(string.ascii_lowercase + string.digits, k=5))
    return f"bedbot-{random_suffix}"

def save_bucket_name(bucket_name):
    """Save bucket name to a temporary file"""
    global bucket_name_file
    if not bucket_name_file:
        bucket_name_file = tempfile.NamedTemporaryFile(mode='w', delete=False, prefix='bedbot_bucket_', suffix='.txt')
        bucket_name_file.write(bucket_name)
        bucket_name_file.close()
        if DEBUG_MODE:
            logger.info(f"Saved bucket name to: {bucket_name_file.name}")

def load_bucket_name():
    """Load bucket name from temporary file if it exists"""
    global bucket_name_file
    # Only load if we have a reference to the current session's bucket file
    if bucket_name_file and hasattr(bucket_name_file, 'name') and os.path.exists(bucket_name_file.name):
        try:
            with open(bucket_name_file.name, 'r') as f:
                bucket_name = f.read().strip()
                if bucket_name:
                    if DEBUG_MODE:
                        logger.info(f"Loaded existing bucket name: {bucket_name}")
                    return bucket_name
        except Exception as e:
            logger.error(f"Error reading bucket name file {bucket_name_file.name}: {e}")
    return None

def cleanup_bucket_name_file():
    """Clean up the bucket name file"""
    global bucket_name_file
    if bucket_name_file and hasattr(bucket_name_file, 'name') and os.path.exists(bucket_name_file.name):
        try:
            os.unlink(bucket_name_file.name)
            if DEBUG_MODE:
                logger.info(f"Cleaned up bucket name file: {bucket_name_file.name}")
        except Exception as e:
            logger.error(f"Error cleaning up bucket name file: {e}")

def create_s3_bucket():
    """Create S3 bucket for file storage"""
    global s3_bucket_name, USE_S3_BUCKET
    if not USE_S3_BUCKET or not s3_client:
        return True
    
    try:
        s3_bucket_name = generate_bucket_name()
        
        # Create bucket
        if profile_region == 'us-east-1':
            s3_client.create_bucket(Bucket=s3_bucket_name)
        else:
            s3_client.create_bucket(
                Bucket=s3_bucket_name,
                CreateBucketConfiguration={'LocationConstraint': profile_region}
            )
        
        # Set bucket policy to private
        s3_client.put_public_access_block(
            Bucket=s3_bucket_name,
            PublicAccessBlockConfiguration={
                'BlockPublicAcls': True,
                'IgnorePublicAcls': True,
                'BlockPublicPolicy': True,
                'RestrictPublicBuckets': True
            }
        )
        
        # Enable server-side encryption with AWS managed keys (AES256)
        s3_client.put_bucket_encryption(
            Bucket=s3_bucket_name,
            ServerSideEncryptionConfiguration={
                'Rules': [
                    {
                        'ApplyServerSideEncryptionByDefault': {
                            'SSEAlgorithm': 'AES256'
                        },
                        'BucketKeyEnabled': True
                    }
                ]
            }
        )
        
        # Save bucket name for Flask restart
        save_bucket_name(s3_bucket_name)
        
        if DEBUG_MODE:
            logger.info(f"Created S3 bucket: {s3_bucket_name}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create S3 bucket: {e}")
        logger.warning("Falling back to local filesystem mode due to S3 bucket creation failure")
        s3_bucket_name = None
        USE_S3_BUCKET = False
        return False

def delete_s3_bucket_with_timeout(timeout_seconds=10):
    """Delete S3 bucket with timeout to prevent hanging on shutdown"""
    def delete_bucket():
        delete_s3_bucket()
    
    thread = threading.Thread(target=delete_bucket)
    thread.daemon = True
    thread.start()
    thread.join(timeout=timeout_seconds)
    
    if thread.is_alive():
        logger.warning(f"S3 bucket cleanup timed out after {timeout_seconds} seconds - bucket may need manual cleanup")

def delete_s3_bucket():
    """Delete S3 bucket and all its contents"""
    global s3_bucket_name
    if not USE_S3_BUCKET or not s3_client or not s3_bucket_name:
        if DEBUG_MODE:
            logger.info("Skipping S3 bucket deletion - not using S3 or no bucket name")
        return
    
    try:
        logger.info(f"Cleaning up S3 bucket: {s3_bucket_name}")
        
        # Check if bucket exists first
        try:
            s3_client.head_bucket(Bucket=s3_bucket_name)
            if DEBUG_MODE:
                logger.info(f"S3 bucket {s3_bucket_name} exists, proceeding with deletion")
        except ClientError as e:
            error_code = e.response['Error']['Code']
            if error_code == '404':
                if DEBUG_MODE:
                    logger.info(f"S3 bucket {s3_bucket_name} doesn't exist - already deleted")
                return
            else:
                raise e
        
        # Delete all objects in bucket with timeout
        paginator = s3_client.get_paginator('list_objects_v2')
        object_count = 0
        
        # Use shorter config for cleanup operations
        cleanup_config = Config(
            read_timeout=5,
            connect_timeout=5,
            retries={'max_attempts': 1}
        )
        cleanup_s3_client = session_aws.client('s3', region_name=profile_region, config=cleanup_config)
        
        for page in paginator.paginate(Bucket=s3_bucket_name):
            if 'Contents' in page:
                objects = [{'Key': obj['Key']} for obj in page['Contents']]
                object_count += len(objects)
                cleanup_s3_client.delete_objects(
                    Bucket=s3_bucket_name,
                    Delete={'Objects': objects}
                )
        
        if object_count > 0:
            logger.info(f"Deleted {object_count} objects from bucket")
        
        # Delete bucket
        cleanup_s3_client.delete_bucket(Bucket=s3_bucket_name)
        logger.info(f"Successfully deleted S3 bucket: {s3_bucket_name}")
        
    except Exception as e:
        logger.warning(f"S3 bucket cleanup failed (bucket may need manual cleanup): {e}")

# Signal handler for graceful shutdown
def signal_handler(signum, frame):
    """Handle interrupt signals for graceful shutdown"""
    signal_name = 'SIGINT' if signum == signal.SIGINT else f'Signal {signum}'
    logger.info(f"Received {signal_name} - shutting down gracefully...")
    cleanup_resources()
    sys.exit(0)

# Register cleanup function to remove temp directories and S3 bucket on shutdown
def cleanup_resources():
    global temp_session_dir, session_upload_folders, cleanup_performed
    
    # Only run cleanup in the main process, not in Flask's debug reloader
    if os.environ.get('WERKZEUG_RUN_MAIN'):
        return
    
    # Prevent duplicate cleanup
    if cleanup_performed:
        return
    cleanup_performed = True
    
    logger.info("Starting cleanup process...")
    
    try:
        # Clean up S3 bucket if using S3 mode (with timeout)
        if USE_S3_BUCKET:
            delete_s3_bucket_with_timeout(timeout_seconds=8)
        
        # Clean up bucket name file
        cleanup_bucket_name_file()
        
        # Clean up all session upload folders (local mode only)
        for folder in list(session_upload_folders):
            try:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
            except Exception as e:
                logger.warning(f"Could not remove session folder {folder}: {e}")
        
        # Clean up Flask session directory
        if temp_session_dir and os.path.exists(temp_session_dir):
            try:
                shutil.rmtree(temp_session_dir)
            except Exception as e:
                logger.warning(f"Could not remove session directory {temp_session_dir}: {e}")
        
        logger.info("Cleanup completed")
        
    except Exception as e:
        logger.error(f"Error during cleanup: {e}")

def initialize_aws_clients():
    """Initialize AWS clients with proper configuration"""
    global bedrock_client, s3_client, session_aws, profile_region
    
    try:
        # Create a session that uses the current AWS profile
        session_aws = boto3.Session()
        
        # Get the region from the profile configuration
        profile_region = session_aws.region_name
        if not profile_region:
            # Fallback to default region if not set in profile
            profile_region = 'us-east-1'
            logger.warning("No region found in AWS profile, defaulting to us-east-1")
        
        # Configure clients with increased timeout for large document processing
        config = Config(
            read_timeout=60,  # 1 minute read timeout (further reduced)
            connect_timeout=30,  # 30 seconds connect timeout
            retries={'max_attempts': 2}  # Reduce retries to 2
        )
        
        bedrock_client = session_aws.client('bedrock-runtime', region_name=profile_region, config=config)
        s3_client = session_aws.client('s3', region_name=profile_region, config=config) if USE_S3_BUCKET else None
        
        if DEBUG_MODE:
            logger.info(f"Initialized Bedrock client with profile: {session_aws.profile_name or 'default'}")
            logger.info(f"Using region: {profile_region}")
            logger.info("Configured clients with 15-minute read timeout for large document processing")
            if USE_S3_BUCKET:
                logger.info("S3 bucket mode enabled")
            else:
                logger.info("Local filesystem mode enabled (--no-bucket)")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to initialize AWS clients: {e}")
        bedrock_client = None
        s3_client = None
        return False

def handle_flask_restart_bucket():
    """Handle S3 bucket name loading during Flask debug restart"""
    global s3_bucket_name, bucket_name_file, USE_S3_BUCKET
    
    # During Flask debug restart, find and load the existing bucket name
    temp_dir = tempfile.gettempdir()
    bucket_files = []
    
    for filename in os.listdir(temp_dir):
        if filename.startswith('bedbot_bucket_') and filename.endswith('.txt'):
            filepath = os.path.join(temp_dir, filename)
            try:
                mtime = os.path.getmtime(filepath)
                bucket_files.append((mtime, filepath))
            except Exception as e:
                logger.error(f"Error checking bucket name file {filepath}: {e}")
    
    if bucket_files:
        # Sort by modification time (newest first) and use the most recent
        bucket_files.sort(reverse=True)
        most_recent_file = bucket_files[0][1]
        
        try:
            with open(most_recent_file, 'r') as f:
                s3_bucket_name = f.read().strip()
                if s3_bucket_name:
                    # Set the bucket_name_file reference for cleanup
                    bucket_name_file = type('obj', (object,), {'name': most_recent_file})()
                    if DEBUG_MODE:
                        logger.info(f"Flask debug restart detected - reusing existing S3 bucket: {s3_bucket_name}")
                else:
                    logger.warning("Flask debug restart detected but bucket name file is empty")
                    USE_S3_BUCKET = False
        except Exception as e:
            logger.error(f"Error reading bucket name file {most_recent_file}: {e}")
            USE_S3_BUCKET = False
    else:
        logger.warning("Flask debug restart detected but no existing bucket name found")
        USE_S3_BUCKET = False

# Allowed file extensions
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'doc', 'docx', 'md', 'json', 'csv'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def truncate_pdf_with_head_and_tail(input_path, output_path, total_mb_size):
    """Truncate a file by taking first 4MB and last 0.5MB"""
    # Constants
    TAIL_SIZE_MB = 0.5
    HEAD_SIZE_MB = total_mb_size - TAIL_SIZE_MB
    
    bytes_head = int(HEAD_SIZE_MB * 1024 * 1024)
    bytes_tail = int(TAIL_SIZE_MB * 1024 * 1024)
    
    try:
        # Get file size
        file_size = os.path.getsize(input_path)
        
        # Handle case where file is smaller than requested size
        if file_size <= (bytes_head + bytes_tail):
            with open(input_path, 'rb') as input_file:
                content = input_file.read()
            
            with open(output_path, 'wb') as output_file:
                output_file.write(content)
            
            logger.info(f"File is smaller than {total_mb_size} MB. Copied entire file.")
            return True
        
        with open(input_path, 'rb') as input_file:
            # Get content from beginning
            head_content = input_file.read(bytes_head)
            
            # Move to position to read tail
            input_file.seek(-bytes_tail, os.SEEK_END)
            tail_content = input_file.read()
        
        # Write combined content
        with open(output_path, 'wb') as output_file:
            output_file.write(head_content)
            output_file.write(tail_content)
        
        actual_size = (len(head_content) + len(tail_content)) / (1024 * 1024)
        logger.info(f"Successfully extracted {HEAD_SIZE_MB:.2f} MB from beginning and {TAIL_SIZE_MB:.2f} MB from end")
        logger.info(f"Output file: {output_path} ({actual_size:.2f} MB)")
        
    except Exception as e:
        logger.error(f"Error truncating file: {e}")
        return False
    
    return True

def merge_pdfs(pdf_paths, output_path, source_directory_name=None):
    """Merge multiple PDF files into a single PDF"""
    try:
        merged_doc = fitz.open()
        
        for pdf_path in pdf_paths:
            try:
                doc = fitz.open(pdf_path)
                merged_doc.insert_pdf(doc)
                doc.close()
                logger.info(f"Added {pdf_path} to merged PDF")
            except Exception as e:
                logger.error(f"Error adding {pdf_path} to merged PDF: {e}")
                continue
        
        merged_doc.save(output_path)
        merged_doc.close()
        
        merged_size = os.path.getsize(output_path) / (1024 * 1024)
        logger.info(f"Successfully merged {len(pdf_paths)} PDFs into {output_path} ({merged_size:.2f} MB)")
        
        return True
        
    except Exception as e:
        logger.error(f"Error merging PDFs: {e}")
        return False

def get_merged_filename_from_files(file_paths):
    """Generate merged filename from first 5 characters of each file, separated by hyphens"""
    if not file_paths:
        return "other"
    
    try:
        # Get all the filenames
        filenames = []
        for file_path in file_paths:
            if hasattr(file_path, 'filename'):
                filenames.append(file_path.filename)
            else:
                filenames.append(str(file_path))
        
        # Extract first 5 characters from each filename (without extension)
        name_parts = []
        for filename in filenames:
            # Remove extension and get base name
            base_name = os.path.splitext(filename)[0]
            # Replace spaces with hyphens and take first 10 characters
            clean_name = base_name.replace(' ', '-')[:10]
            if clean_name and clean_name not in name_parts:
                name_parts.append(clean_name)
        
        # Join with hyphens
        if name_parts:
            merged_name = '-'.join(name_parts)
            # Ensure no double hyphens and clean up
            merged_name = '-'.join(part for part in merged_name.split('-') if part)
            return merged_name if merged_name else "other"
        
        # Fallback to "other"
        return "other"
        
    except Exception:
        return "other"

def convert_pdf_to_markdown_bedrock(pdf_path_or_s3key, is_s3=False, s3_bucket=None, model_id=None):
    """Convert a PDF file to markdown using Bedrock model"""
    if not bedrock_client:
        logger.error("Bedrock client not initialized")
        return None
    
    if not model_id:
        model_id = PDF_CONVERT_MODEL
    
    # Retry logic for rate limiting
    max_retries = 3
    base_delay = 2
    
    for attempt in range(max_retries):
        try:
            # Prepare the prompt for markdown conversion
            conversion_prompt = """Please convert this PDF document to markdown format. 

IMPORTANT INSTRUCTIONS:
- Extract all text content from the document
- Preserve the document structure using appropriate markdown headers (# ## ###)
- Convert tables to markdown table format if present
- Include all relevant text content
- Format lists using markdown list syntax
- Preserve paragraph breaks and spacing
- If there are multiple sections, use appropriate header levels
- Do not add any commentary or explanation - just provide the markdown conversion

Return only the markdown content without any additional text."""
            
            content = []
            
            if is_s3:
                # S3 mode - try S3 location first, fall back to bytes if not supported
                if not s3_bucket:
                    logger.error("S3 bucket name required for S3 mode")
                    return None
                
                # Sanitize document name for Bedrock
                doc_name = os.path.basename(pdf_path_or_s3key).replace('.pdf', '')
                doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
                doc_name = re.sub(r'\s+', ' ', doc_name).strip()
                doc_name = doc_name[:50]
                
                # Try S3 location first (for models that support it)
                try:
                    # Get AWS account ID for bucket owner
                    account_id = session_aws.client('sts').get_caller_identity()['Account']
                    
                    logger.info(f"Converting PDF to markdown via S3 location: {doc_name} (s3://{s3_bucket}/{pdf_path_or_s3key})")
                    
                    content.append({
                        "document": {
                            "format": "pdf",
                            "name": doc_name,
                            "source": {
                                "s3Location": {
                                    "uri": f"s3://{s3_bucket}/{pdf_path_or_s3key}",
                                    "bucketOwner": account_id
                                }
                            }
                        }
                    })
                    
                    # Try the conversion with S3 location
                    content.append({
                        "text": conversion_prompt
                    })
                    
                    messages = [{
                        "role": "user",
                        "content": content
                    }]
                    
                    inference_config = {
                        "maxTokens": 4000,
                        "temperature": 0.3,
                        "topP": 0.9
                    }
                    
                    response = bedrock_client.converse(
                        modelId=model_id,
                        messages=messages,
                        inferenceConfig=inference_config
                    )
                    
                    markdown_content = response['output']['message']['content'][0]['text']
                    logger.info(f"Successfully converted PDF to markdown using S3 location ({len(markdown_content)} characters)")
                    return markdown_content
                    
                except ClientError as e:
                    error_code = e.response.get('Error', {}).get('Code', 'Unknown')
                    if 'ValidationException' in str(e) and 'S3Location' in str(e):
                        logger.info(f"Model {model_id} doesn't support S3Location, falling back to bytes method")
                        # Fall back to downloading and sending as bytes
                        try:
                            response = s3_client.get_object(Bucket=s3_bucket, Key=pdf_path_or_s3key)
                            pdf_bytes = response['Body'].read()
                            
                            logger.info(f"Converting PDF to markdown (downloaded from S3): {doc_name} ({len(pdf_bytes)} bytes)")
                            
                            content = [{
                                "document": {
                                    "format": "pdf",
                                    "name": doc_name,
                                    "source": {
                                        "bytes": pdf_bytes
                                    }
                                }
                            }]
                        except Exception as download_error:
                            logger.error(f"Failed to download PDF from S3: {download_error}")
                            return None
                    else:
                        # Other error, re-raise for retry logic
                        raise e
            else:
                # Local mode - read file bytes
                if not os.path.exists(pdf_path_or_s3key):
                    logger.error(f"PDF file not found: {pdf_path_or_s3key}")
                    return None
                    
                with open(pdf_path_or_s3key, 'rb') as f:
                    pdf_bytes = f.read()
                
                # Sanitize document name for Bedrock
                doc_name = os.path.basename(pdf_path_or_s3key).replace('.pdf', '')
                doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
                doc_name = re.sub(r'\s+', ' ', doc_name).strip()
                doc_name = doc_name[:50]
                
                logger.info(f"Converting PDF to markdown: {doc_name} ({len(pdf_bytes)} bytes)")
                
                content.append({
                    "document": {
                        "format": "pdf",
                        "name": doc_name,
                        "source": {
                            "bytes": pdf_bytes
                        }
                    }
                })
            
            content.append({
                "text": conversion_prompt
            })
            
            messages = [{
                "role": "user",
                "content": content
            }]
            
            inference_config = {
                "maxTokens": 4000,
                "temperature": 0.3,  # Lower temperature for more consistent conversion
                "topP": 0.9
            }
            
            logger.info(f"Sending PDF conversion request to Bedrock model: {model_id} (attempt {attempt + 1})")
            
            response = bedrock_client.converse(
                modelId=model_id,
                messages=messages,
                inferenceConfig=inference_config
            )
            
            markdown_content = response['output']['message']['content'][0]['text']
            logger.info(f"Successfully converted PDF to markdown ({len(markdown_content)} characters)")
            
            return markdown_content
            
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', 'Unknown')
            if error_code in ['ServiceUnavailableException', 'ThrottlingException'] and attempt < max_retries - 1:
                # Rate limiting or service unavailable - wait and retry
                delay = base_delay * (2 ** attempt)  # Exponential backoff
                logger.warning(f"Rate limit hit on attempt {attempt + 1}, retrying in {delay} seconds...")
                import time
                time.sleep(delay)
                continue
            else:
                logger.error(f"Bedrock API error during PDF conversion: {e}")
                return None
        except Exception as e:
            logger.error(f"Unexpected error during PDF conversion: {e}")
            return None
    
    logger.error(f"Failed to convert PDF after {max_retries} attempts")
    return None

def convert_pdfs_parallel(pdf_conversion_tasks):
    """Convert multiple PDFs to markdown in parallel using ThreadPoolExecutor"""
    if not pdf_conversion_tasks:
        return {}
    
    logger.info(f"Starting parallel PDF conversion for {len(pdf_conversion_tasks)} files (max concurrent: {BEDROCK_MAX_CONCURRENT})")
    
    results = {}
    
    def convert_single_pdf(task):
        """Convert a single PDF - wrapper for thread pool"""
        file_key, s3_key, s3_bucket, filename = task
        try:
            logger.info(f"[Parallel] Starting conversion: {filename}")
            markdown_content = convert_pdf_to_markdown_bedrock(s3_key, is_s3=True, s3_bucket=s3_bucket, model_id=PDF_CONVERT_MODEL)
            logger.info(f"[Parallel] Completed conversion: {filename} ({len(markdown_content) if markdown_content else 0} chars)")
            return file_key, markdown_content, None
        except Exception as e:
            logger.error(f"[Parallel] Failed conversion: {filename} - {e}")
            return file_key, None, str(e)
    
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=BEDROCK_MAX_CONCURRENT) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(convert_single_pdf, task): task 
            for task in pdf_conversion_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_task):
            task = future_to_task[future]
            file_key, markdown_content, error = future.result()
            
            results[file_key] = {
                'markdown_content': markdown_content,
                'error': error,
                'filename': task[3]  # filename is the 4th element
            }
    
    logger.info(f"Completed parallel PDF conversion. Successful: {sum(1 for r in results.values() if r['markdown_content'])}, Failed: {sum(1 for r in results.values() if r['error'])}")
    return results

def read_file_content(filepath_or_s3key, is_s3=False):
    """Read content from uploaded file or return file info for PDF processing"""
    try:
        if is_s3:
            # S3 mode - filepath_or_s3key is the S3 key
            filename = os.path.basename(filepath_or_s3key).lower()
            
            if filename.endswith('.pdf'):
                # For PDFs, return file info for Bedrock document processing
                return {
                    'type': 'pdf',
                    's3_key': filepath_or_s3key,
                    'filename': os.path.basename(filepath_or_s3key)
                }
            else:
                # For text files, read content from S3
                response = s3_client.get_object(Bucket=s3_bucket_name, Key=filepath_or_s3key)
                return response['Body'].read().decode('utf-8', errors='ignore')
        else:
            # Local mode - filepath_or_s3key is the local file path
            filename = os.path.basename(filepath_or_s3key).lower()
            
            if filename.endswith('.pdf'):
                # For PDFs, return file info for Bedrock document processing
                return {
                    'type': 'pdf',
                    'path': filepath_or_s3key,
                    'filename': os.path.basename(filepath_or_s3key)
                }
            else:
                # For text files, read content normally
                with open(filepath_or_s3key, 'r', encoding='utf-8', errors='ignore') as f:
                    return f.read()
    except Exception as e:
        logger.error(f"Error reading file {filepath_or_s3key}: {e}")
        return None

def call_bedrock(prompt, context="", pdf_files=None):
    """Call AWS Bedrock model using Converse API with document support"""
    if not bedrock_client:
        return "Error: Bedrock client not initialized. Please check your AWS credentials."
    
    try:
        content = []
        
        # Debug logging
        logger.info(f"PDF files provided: {len(pdf_files) if pdf_files else 0}")
        if pdf_files:
            for i, pdf_info in enumerate(pdf_files):
                if USE_S3_BUCKET:
                    logger.info(f"PDF {i+1}: {pdf_info.get('filename', 'unknown')} at s3://{s3_bucket_name}/{pdf_info.get('s3_key', 'unknown')}")
                else:
                    logger.info(f"PDF {i+1}: {pdf_info.get('filename', 'unknown')} at {pdf_info.get('path', 'unknown')}")
        
        # Add PDF documents if provided (prefer markdown when available)
        if pdf_files:
            for pdf_info in pdf_files:
                try:
                    # Check if we have markdown version available
                    has_markdown = pdf_info.get('has_markdown', False)
                    doc_name = pdf_info['filename'].replace('.pdf', '')
                    doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
                    doc_name = re.sub(r'\s+', ' ', doc_name).strip()
                    doc_name = doc_name[:50]
                    
                    if has_markdown:
                        # Use markdown content instead of PDF for faster processing
                        if USE_S3_BUCKET and 'md_s3_key' in pdf_info:
                            # S3 mode - read markdown from S3
                            try:
                                response = s3_client.get_object(Bucket=s3_bucket_name, Key=pdf_info['md_s3_key'])
                                markdown_content = response['Body'].read().decode('utf-8', errors='ignore')
                                logger.info(f"Using markdown content from S3 for faster processing: {doc_name} (s3://{s3_bucket_name}/{pdf_info['md_s3_key']})")
                                
                                # Add markdown content as text instead of document
                                content.append({
                                    "text": f"\n\n--- Content from {pdf_info['filename']} (PDF converted to markdown) ---\n{markdown_content}"
                                })
                                continue
                            except Exception as e:
                                logger.warning(f"Failed to read markdown from S3, falling back to PDF: {e}")
                                has_markdown = False
                        elif 'md_path' in pdf_info and os.path.exists(pdf_info['md_path']):
                            # Local mode - read markdown from file
                            try:
                                with open(pdf_info['md_path'], 'r', encoding='utf-8', errors='ignore') as md_file:
                                    markdown_content = md_file.read()
                                logger.info(f"Using markdown content from file for faster processing: {doc_name} ({pdf_info['md_path']})")
                                
                                # Add markdown content as text instead of document
                                content.append({
                                    "text": f"\n\n--- Content from {pdf_info['filename']} (PDF converted to markdown) ---\n{markdown_content}"
                                })
                                continue
                            except Exception as e:
                                logger.warning(f"Failed to read markdown from file, falling back to PDF: {e}")
                                has_markdown = False
                    
                    # Fall back to PDF processing if no markdown or markdown read failed
                    if USE_S3_BUCKET and 's3_key' in pdf_info:
                        # S3 mode - try S3 location first, fall back to bytes if model doesn't support it
                        s3_key = pdf_info['s3_key']
                        
                        # Check if this model supports S3Location by trying a simple test
                        # For efficiency, we'll just try to download and use bytes for non-Nova models
                        if 'nova' in BEDROCK_MODEL.lower():
                            # Nova models support S3Location
                            try:
                                # Get AWS account ID for bucket owner
                                account_id = session_aws.client('sts').get_caller_identity()['Account']
                                
                                logger.info(f"Adding PDF document from S3 location: {doc_name} (s3://{s3_bucket_name}/{s3_key})")
                                
                                content.append({
                                    "document": {
                                        "format": "pdf",
                                        "name": doc_name,
                                        "source": {
                                            "s3Location": {
                                                "uri": f"s3://{s3_bucket_name}/{s3_key}",
                                                "bucketOwner": account_id
                                            }
                                        }
                                    }
                                })
                            except Exception as e:
                                logger.warning(f"Failed to get account ID, falling back to bytes method: {e}")
                                # Fall back to bytes method
                                try:
                                    response = s3_client.get_object(Bucket=s3_bucket_name, Key=s3_key)
                                    pdf_bytes = response['Body'].read()
                                    
                                    logger.info(f"Adding PDF document (downloaded from S3): {doc_name} ({len(pdf_bytes)} bytes)")
                                    
                                    content.append({
                                        "document": {
                                            "format": "pdf",
                                            "name": doc_name,
                                            "source": {
                                                "bytes": pdf_bytes
                                            }
                                        }
                                    })
                                except Exception as download_error:
                                    logger.error(f"Failed to download PDF from S3: {download_error}")
                                    continue
                        else:
                            # Non-Nova models - use bytes method directly
                            try:
                                response = s3_client.get_object(Bucket=s3_bucket_name, Key=s3_key)
                                pdf_bytes = response['Body'].read()
                                
                                logger.info(f"Adding PDF document (downloaded from S3): {doc_name} ({len(pdf_bytes)} bytes)")
                                
                                content.append({
                                    "document": {
                                        "format": "pdf",
                                        "name": doc_name,
                                        "source": {
                                            "bytes": pdf_bytes
                                        }
                                    }
                                })
                            except Exception as download_error:
                                logger.error(f"Failed to download PDF from S3: {download_error}")
                                continue
                    else:
                        # Local mode - read file bytes
                        pdf_path = pdf_info['path']
                        if not os.path.exists(pdf_path):
                            logger.error(f"PDF file not found: {pdf_path}")
                            continue
                            
                        with open(pdf_path, 'rb') as f:
                            pdf_bytes = f.read()
                        
                        logger.info(f"Adding PDF document: {doc_name} ({len(pdf_bytes)} bytes)")
                        
                        content.append({
                            "document": {
                                "format": "pdf",
                                "name": doc_name,
                                "source": {
                                    "bytes": pdf_bytes
                                }
                            }
                        })
                except Exception as e:
                    logger.error(f"Error processing PDF {pdf_info}: {e}")
        
        # Add text context if available and no PDFs
        if context and not pdf_files:
            enhanced_prompt = f"{context}\n\n{prompt}"
        else:
            # Create enhanced prompt for document analysis
            # Count total original documents uploaded by user
            total_uploaded_docs = len(session.get('uploaded_files', []))
            
            if pdf_files and len(pdf_files) > 1:
                # When multiple PDFs are provided, give clear instructions for comparison
                enhanced_prompt = textwrap.dedent(f"""I have provided {total_uploaded_docs} documents for analysis. Please carefully analyze ALL the provided documents and answer the following question:

                {prompt}

                IMPORTANT INSTRUCTIONS:
                - You have access to {total_uploaded_docs} documents - please read and analyze ALL of them
                - Compare and analyze the content across all provided documents
                - If one document contains requirements or criteria, evaluate the other documents against those criteria
                - Provide specific names, examples, and evidence from the documents to support your analysis
                - Be thorough and analytical in your comparison
                - When asked for the option or to compare options presented in different documents, provide specific and detailed comparisons""")
            elif pdf_files and len(pdf_files) == 1:
                enhanced_prompt = textwrap.dedent(f"""I have provided {total_uploaded_docs} document(s) for analysis. Please carefully analyze the provided document(s) and answer the following question:

                {prompt}

                Please provide specific information and examples from the document(s) to support your response.""")
            else:
                enhanced_prompt = prompt
            
        logger.info(f"Enhanced prompt length: {len(enhanced_prompt)} characters")
        logger.info(f"Total content items: {len(content) + 1} (PDFs: {len(content)}, text: 1)")
            
        content.append({
            "text": enhanced_prompt
        })
        
        messages = [
            {
                "role": "user",
                "content": content
            }
        ]
        
        inference_config = {
            "maxTokens": 4000,
            "temperature": 0.7,
            "topP": 0.9
        }
        
        logger.info(f"Sending request to Bedrock with {len(content)} content items")
        
        # Debug: Print messages JSON structure if debug mode is enabled
        if DEBUG_MODE:
            logger.info("=== DEBUG: Converse API Messages Structure ===")
            logger.info(json.dumps(messages, indent=2, default=str))
            logger.info("=== END DEBUG ===")
        
        try:
            response = bedrock_client.converse(
                modelId=BEDROCK_MODEL,
                messages=messages,
                inferenceConfig=inference_config
            )
            
            return response['output']['message']['content'][0]['text']
            
        except ClientError as api_error:
            # Check if this is an S3Location validation error
            if ('ValidationException' in str(api_error) and 'S3Location' in str(api_error) and 
                USE_S3_BUCKET and pdf_files):
                logger.warning(f"Model {BEDROCK_MODEL} doesn't support S3Location, retrying with bytes method")
                
                # Rebuild content with bytes method for PDFs
                retry_content = []
                
                # Re-add PDF documents using bytes method
                for pdf_info in pdf_files:
                    try:
                        # Check if we have markdown version available (prefer markdown)
                        has_markdown = pdf_info.get('has_markdown', False)
                        doc_name = pdf_info['filename'].replace('.pdf', '')
                        doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
                        doc_name = re.sub(r'\s+', ' ', doc_name).strip()
                        doc_name = doc_name[:50]
                        
                        if has_markdown:
                            # Use markdown content
                            if USE_S3_BUCKET and 'md_s3_key' in pdf_info:
                                try:
                                    response = s3_client.get_object(Bucket=s3_bucket_name, Key=pdf_info['md_s3_key'])
                                    markdown_content = response['Body'].read().decode('utf-8', errors='ignore')
                                    retry_content.append({
                                        "text": f"\n\n--- Content from {pdf_info['filename']} (PDF converted to markdown) ---\n{markdown_content}"
                                    })
                                    continue
                                except Exception as md_error:
                                    logger.warning(f"Failed to read markdown in retry, falling back to PDF bytes: {md_error}")
                            elif 'md_path' in pdf_info and os.path.exists(pdf_info['md_path']):
                                try:
                                    with open(pdf_info['md_path'], 'r', encoding='utf-8', errors='ignore') as md_file:
                                        markdown_content = md_file.read()
                                    retry_content.append({
                                        "text": f"\n\n--- Content from {pdf_info['filename']} (PDF converted to markdown) ---\n{markdown_content}"
                                    })
                                    continue
                                except Exception as md_error:
                                    logger.warning(f"Failed to read markdown file in retry, falling back to PDF bytes: {md_error}")
                        
                        # Fall back to PDF bytes method
                        if USE_S3_BUCKET and 's3_key' in pdf_info:
                            try:
                                response = s3_client.get_object(Bucket=s3_bucket_name, Key=pdf_info['s3_key'])
                                pdf_bytes = response['Body'].read()
                                
                                logger.info(f"Retry: Adding PDF document (downloaded from S3): {doc_name} ({len(pdf_bytes)} bytes)")
                                
                                retry_content.append({
                                    "document": {
                                        "format": "pdf",
                                        "name": doc_name,
                                        "source": {
                                            "bytes": pdf_bytes
                                        }
                                    }
                                })
                            except Exception as download_error:
                                logger.error(f"Retry: Failed to download PDF from S3: {download_error}")
                                continue
                        else:
                            # Local mode
                            pdf_path = pdf_info['path']
                            if os.path.exists(pdf_path):
                                with open(pdf_path, 'rb') as f:
                                    pdf_bytes = f.read()
                                
                                retry_content.append({
                                    "document": {
                                        "format": "pdf",
                                        "name": doc_name,
                                        "source": {
                                            "bytes": pdf_bytes
                                        }
                                    }
                                })
                    except Exception as pdf_error:
                        logger.error(f"Retry: Error processing PDF {pdf_info}: {pdf_error}")
                
                # Add the text prompt
                retry_content.append({
                    "text": enhanced_prompt
                })
                
                # Retry the API call with bytes method
                retry_messages = [{
                    "role": "user",
                    "content": retry_content
                }]
                
                logger.info(f"Retrying request with {len(retry_content)} content items using bytes method")
                
                retry_response = bedrock_client.converse(
                    modelId=BEDROCK_MODEL,
                    messages=retry_messages,
                    inferenceConfig=inference_config
                )
                
                return retry_response['output']['message']['content'][0]['text']
            else:
                # Other API error, re-raise
                raise api_error
        
    except ClientError as e:
        logger.error(f"Bedrock API error: {e}")
        return f"Error calling Bedrock: {str(e)}"
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return f"Unexpected error: {str(e)}"

def convert_markdown_to_html(text):
    """Convert markdown text to HTML with enhanced formatting"""
    # Configure markdown with useful extensions
    md = markdown.Markdown(extensions=[
        'extra',          # Includes tables, fenced code blocks, etc.
        'codehilite',     # Syntax highlighting for code blocks
        'toc',            # Table of contents
        'nl2br',          # Convert newlines to <br> tags
        'sane_lists'      # Better list handling
    ])
    
    return md.convert(text)

def get_session_upload_location():
    """Get or create session-specific upload location (S3 prefix or local folder)"""
    global temp_session_dir, session_upload_folders
    
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    
    # Create Flask session directory if it doesn't exist
    if not temp_session_dir:
        temp_session_dir = tempfile.mkdtemp(prefix='bedbot_sessions_')
        os.chmod(temp_session_dir, 0o700)  # Only owner can read/write/execute
        app.config['SESSION_FILE_DIR'] = temp_session_dir
        if DEBUG_MODE:
            logger.info(f"Created temporary session directory: {temp_session_dir}")
    
    if USE_S3_BUCKET:
        # S3 mode - return S3 prefix for this session
        s3_prefix = f"session_{session['session_id']}/"
        session['s3_prefix'] = s3_prefix
        session.modified = True
        if DEBUG_MODE:
            logger.info(f"Using S3 session prefix: {s3_prefix}")
        return s3_prefix
    else:
        # Local mode - create temporary directory
        session_folder = tempfile.mkdtemp(prefix=f'bedbot_session_{session["session_id"]}_')
        os.chmod(session_folder, 0o700)  # Only owner can read/write/execute
        
        # Track this folder for cleanup on shutdown
        session_upload_folders.add(session_folder)
        
        # Store the session folder path in the session for cleanup
        session['session_folder'] = session_folder
        session.modified = True
        
        if DEBUG_MODE:
            logger.info(f"Created session upload folder: {session_folder}")
        return session_folder

@app.route('/')
def index():
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
        session['chat_history'] = []
        session['document_context'] = ""
        session['uploaded_files'] = []
        session['pdf_files'] = []
        # Create session upload location
        get_session_upload_location()
    return render_template('index.html')

@app.route('/chat', methods=['POST'])
def chat():
    try:
        data = request.get_json()
        user_message = data.get('message', '').strip()
        
        if not user_message:
            return jsonify({'error': 'Empty message'}), 400
        
        # Get document context and PDF files
        context = session.get('document_context', '')
        pdf_files = session.get('pdf_files', [])
        
        # Debug logging
        logger.info(f"Chat request - PDF files in session: {len(pdf_files)}")
        logger.info(f"Chat request - Context length: {len(context)}")
        
        # Always call with PDF files (empty list if none)
        bot_response = call_bedrock(user_message, context, pdf_files)
        
        # Convert markdown to HTML
        bot_response_html = convert_markdown_to_html(bot_response)
        
        # Store in chat history (store both markdown and HTML)
        if 'chat_history' not in session:
            session['chat_history'] = []
        
        session['chat_history'].append({
            'user': user_message,
            'bot': bot_response,  # Store original markdown
            'bot_html': bot_response_html,  # Store HTML version
            'timestamp': datetime.now().isoformat()
        })
        
        session.modified = True
        
        return jsonify({
            'response': bot_response_html,  # Send HTML to frontend
            'timestamp': datetime.now().isoformat()
        })
        
    except Exception as e:
        logger.error(f"Chat error: {e}")
        return jsonify({'error': 'Internal server error'}), 500

# Removed complex session locking

@app.route('/upload', methods=['POST'])
def upload_files():
    try:
        import time
        request_start_time = time.time()
        logger.info(f"=== UPLOAD REQUEST START (timestamp: {request_start_time}) ===")
        if 'files' not in request.files:
            return jsonify({'error': 'No files provided'}), 400
        
        files = request.files.getlist('files')
        
        # Check file count limit
        current_file_count = len(session.get('uploaded_files', []))
        if current_file_count + len(files) > MAX_FILES_PER_SESSION:
            return jsonify({'error': f'Maximum {MAX_FILES_PER_SESSION} files allowed per session'}), 400
        
        uploaded_files = []
        new_text_content = ""
        new_pdf_files = []
        
        # Get or create session-specific upload location
        if USE_S3_BUCKET:
            if 's3_prefix' in session:
                session_location = session['s3_prefix']
            else:
                session_location = get_session_upload_location()
        else:
            if 'session_folder' in session and os.path.exists(session['session_folder']):
                session_location = session['session_folder']
            else:
                session_location = get_session_upload_location()
        
        # Check if we should merge PDFs (S3 mode, multiple files, multiple PDFs)
        pdf_files_to_merge = []
        non_pdf_files = []
        
        for file in files:
            if file and file.filename and allowed_file(file.filename):
                if file.filename.lower().endswith('.pdf'):
                    pdf_files_to_merge.append(file)
                else:
                    non_pdf_files.append(file)
        
        should_merge_pdfs = (PDF_MERGE_ENABLED and 
                           USE_S3_BUCKET and 
                           len(pdf_files_to_merge) > 1 and 
                           len(files) > 1)
        
        if DEBUG_MODE and len(pdf_files_to_merge) > 1:
            if not PDF_MERGE_ENABLED:
                logger.info(f"PDF merging disabled via PDF_MERGE=0 - processing {len(pdf_files_to_merge)} PDFs individually")
            elif not USE_S3_BUCKET:
                logger.info("PDF merging only available in S3 mode - processing PDFs individually")
        
        if should_merge_pdfs:
            logger.info(f"Merging {len(pdf_files_to_merge)} PDF files for S3 upload")
            
            # Save PDFs temporarily for merging
            temp_pdf_paths = []
            for pdf_file in pdf_files_to_merge:
                temp_pdf = tempfile.NamedTemporaryFile(delete=False, suffix='.pdf')
                pdf_file.save(temp_pdf.name)
                temp_pdf.close()
                temp_pdf_paths.append(temp_pdf.name)
            
            # Generate merged filename
            merged_name = get_merged_filename_from_files(pdf_files_to_merge)
            random_suffix = ''.join(random.choices(string.ascii_lowercase + string.digits, k=5))
            merged_filename = f"{merged_name}-{random_suffix}.pdf"
            timestamped_merged_filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{merged_filename}"
            
            # Create merged PDF
            temp_merged = tempfile.NamedTemporaryFile(delete=False, suffix='.pdf')
            temp_merged.close()
            
            if merge_pdfs(temp_pdf_paths, temp_merged.name, merged_name):
                # In S3 mode, don't truncate - allow larger files
                merged_size = os.path.getsize(temp_merged.name)
                logger.info(f'Merged PDF is {merged_size / (1024*1024):.1f} MB (no truncation in S3 mode)')
                
                # Upload merged PDF to S3
                s3_key = f"{session_location}{timestamped_merged_filename}"
                
                try:
                    with open(temp_merged.name, 'rb') as merged_file:
                        s3_client.upload_fileobj(merged_file, s3_bucket_name, s3_key)
                    logger.info(f"Uploaded merged PDF to S3: s3://{s3_bucket_name}/{s3_key}")
                    
                    # Convert merged PDF to markdown
                    logger.info(f"Converting merged PDF to markdown: {merged_filename}")
                    markdown_content = convert_pdf_to_markdown_bedrock(s3_key, is_s3=True, s3_bucket=s3_bucket_name, model_id=PDF_CONVERT_MODEL)
                    
                    if markdown_content:
                        # Create markdown filename
                        base_name = os.path.splitext(timestamped_merged_filename)[0]
                        md_filename = f"{base_name}.md"
                        md_s3_key = f"{session_location}{md_filename}"
                        
                        # Upload markdown to S3
                        try:
                            s3_client.put_object(
                                Bucket=s3_bucket_name,
                                Key=md_s3_key,
                                Body=markdown_content.encode('utf-8'),
                                ContentType='text/markdown'
                            )
                            logger.info(f"Uploaded merged PDF markdown to S3: s3://{s3_bucket_name}/{md_s3_key}")
                            
                            # Add merged PDF with markdown to results
                            new_pdf_files.append({
                                's3_key': s3_key,
                                'md_s3_key': md_s3_key,
                                'filename': merged_filename,
                                'timestamped_filename': timestamped_merged_filename,
                                'has_markdown': True
                            })
                            
                            # Also add markdown content to text context for immediate use
                            new_text_content += f"\n\n--- Content from {merged_filename} (merged PDF converted to markdown) ---\n{markdown_content}"
                            
                        except Exception as e:
                            logger.error(f"Error uploading merged PDF markdown to S3: {e}")
                            # Fall back to PDF-only storage
                            new_pdf_files.append({
                                's3_key': s3_key,
                                'filename': merged_filename,
                                'timestamped_filename': timestamped_merged_filename,
                                'has_markdown': False
                            })
                    else:
                        logger.warning(f"Failed to convert merged PDF to markdown: {merged_filename}")
                        # Add merged PDF without markdown
                        new_pdf_files.append({
                            's3_key': s3_key,
                            'filename': merged_filename,
                            'timestamped_filename': timestamped_merged_filename,
                            'has_markdown': False
                        })
                    
                    uploaded_files.append({
                        'filename': merged_filename,
                        'size': merged_size,
                        'status': 'success',
                        'type': 'pdf'
                    })
                    
                except Exception as e:
                    logger.error(f"Error uploading merged PDF to S3: {e}")
                    uploaded_files.append({
                        'filename': merged_filename,
                        'status': 'error',
                        'message': 'Failed to upload merged PDF to S3'
                    })
            
            # Clean up temporary files
            for temp_path in temp_pdf_paths:
                os.unlink(temp_path)
            os.unlink(temp_merged.name)
            
            # Process non-PDF files normally
            files_to_process = non_pdf_files
        else:
            # Process all files normally
            files_to_process = files
        
        logger.info(f"Processing {len(files_to_process)} files individually")
        
        # Collect PDFs for parallel processing
        pdf_conversion_tasks = []
        pdf_files_info = {}  # Track PDF info by file key
        
        for file in files_to_process:
            if file and file.filename and allowed_file(file.filename):
                logger.info(f"Processing file: {file.filename}")
                # Check file size
                file.seek(0, 2)  # Seek to end
                file_size = file.tell()
                file.seek(0)  # Reset to beginning
                
                needs_truncation = file_size > MAX_FILE_SIZE and not USE_S3_BUCKET
                if needs_truncation:
                    logger.info(f'File "{file.filename}" is {file_size / (1024*1024):.1f} MB, truncating to {MAX_FILE_SIZE / (1024*1024):.1f} MB')
                elif USE_S3_BUCKET and file_size > MAX_FILE_SIZE:
                    logger.info(f'File "{file.filename}" is {file_size / (1024*1024):.1f} MB (no truncation in S3 mode)')
                
                filename = secure_filename(file.filename)
                # Add timestamp to avoid conflicts
                timestamped_filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{filename}"
                
                if USE_S3_BUCKET:
                    # S3 mode - upload to S3
                    s3_key = f"{session_location}{timestamped_filename}"
                    
                    try:
                        if needs_truncation:
                            # Save file temporarily for truncation
                            temp_input = tempfile.NamedTemporaryFile(delete=False)
                            file.save(temp_input.name)
                            temp_input.close()
                            
                            # Create truncated version
                            temp_output = tempfile.NamedTemporaryFile(delete=False)
                            temp_output.close()
                            
                            if truncate_pdf_with_head_and_tail(temp_input.name, temp_output.name, MAX_FILE_SIZE / (1024*1024)):
                                # Upload truncated file
                                with open(temp_output.name, 'rb') as truncated_file:
                                    s3_client.upload_fileobj(truncated_file, s3_bucket_name, s3_key)
                                logger.info(f"Uploaded truncated file to S3: s3://{s3_bucket_name}/{s3_key}")
                                file_size = os.path.getsize(temp_output.name)  # Update file size
                            else:
                                raise Exception("File truncation failed")
                            
                            # Clean up temp files
                            os.unlink(temp_input.name)
                            os.unlink(temp_output.name)
                        else:
                            # Upload original file
                            s3_client.upload_fileobj(file, s3_bucket_name, s3_key)
                            logger.info(f"Uploaded file to S3: s3://{s3_bucket_name}/{s3_key}")
                        
                        # Handle file content based on type
                        try:
                            content = read_file_content(s3_key, is_s3=True)
                            if content:
                                if isinstance(content, dict) and content.get('type') == 'pdf':
                                    # Always add PDF to the list first, regardless of markdown conversion success
                                    pdf_info = {
                                        's3_key': s3_key,
                                        'filename': file.filename,
                                        'timestamped_filename': timestamped_filename,
                                        'has_markdown': False
                                    }
                                    
                                    # Add to parallel conversion queue if enabled
                                    if PDF_CONVERSION_ENABLED:
                                        file_key = f"{file.filename}_{len(pdf_conversion_tasks)}"  # Unique key
                                        pdf_conversion_tasks.append((file_key, s3_key, s3_bucket_name, file.filename))
                                        pdf_files_info[file_key] = {
                                            'pdf_info': pdf_info,
                                            'session_location': session_location,
                                            'timestamped_filename': timestamped_filename
                                        }
                                        logger.info(f"Added PDF to parallel conversion queue: {file.filename}")
                                    else:
                                        logger.info(f"PDF conversion disabled - PDF will be processed directly by Bedrock: {file.filename}")
                                    
                                    # Always add the PDF file (markdown will be added later if conversion succeeds)
                                    new_pdf_files.append(pdf_info)
                                    logger.info(f"Added PDF to new_pdf_files: {pdf_info['filename']} (will attempt conversion: {PDF_CONVERSION_ENABLED})")
                                    
                                    uploaded_files.append({
                                        'filename': file.filename,
                                        'size': file_size,
                                        'status': 'success',
                                        'type': 'pdf'
                                    })
                                else:
                                    # Regular text content
                                    new_text_content += f"\n\n--- Content from {file.filename} ---\n{content}"
                                    uploaded_files.append({
                                        'filename': file.filename,
                                        'size': file_size,
                                        'status': 'success',
                                        'type': 'text'
                                    })
                            else:
                                uploaded_files.append({
                                    'filename': file.filename,
                                    'status': 'error',
                                    'message': 'Could not read file content'
                                })
                        except Exception as e:
                            logger.error(f"Error processing file content for {file.filename}: {e}")
                            import traceback
                            logger.error(f"Full traceback: {traceback.format_exc()}")
                            uploaded_files.append({
                                'filename': file.filename,
                                'status': 'error',
                                'message': f'Error processing file: {str(e)}'
                            })
                    except Exception as e:
                        logger.error(f"Error uploading {file.filename} to S3: {e}")
                        import traceback
                        logger.error(f"Full traceback: {traceback.format_exc()}")
                        uploaded_files.append({
                            'filename': file.filename,
                            'status': 'error',
                            'message': 'Failed to upload to S3'
                        })
                    
                    logger.info(f"Completed processing file: {file.filename}")
                else:
                    # Local mode - save to local filesystem
                    filepath = os.path.join(session_location, timestamped_filename)
                    
                    if needs_truncation:
                        # Save file temporarily for truncation
                        temp_input = tempfile.NamedTemporaryFile(delete=False)
                        file.save(temp_input.name)
                        temp_input.close()
                        
                        # Create truncated version directly to final path
                        if truncate_pdf_with_head_and_tail(temp_input.name, filepath, MAX_FILE_SIZE / (1024*1024)):
                            logger.info(f"Saved truncated file to: {filepath}")
                            file_size = os.path.getsize(filepath)  # Update file size
                        else:
                            raise Exception("File truncation failed")
                        
                        # Clean up temp file
                        os.unlink(temp_input.name)
                    else:
                        # Save original file
                        file.save(filepath)
                    
                    # Handle file content based on type
                    content = read_file_content(filepath, is_s3=False)
                    if content:
                        if isinstance(content, dict) and content.get('type') == 'pdf':
                            # Convert PDF to markdown using Bedrock
                            logger.info(f"Converting PDF to markdown: {file.filename}")
                            markdown_content = convert_pdf_to_markdown_bedrock(filepath, is_s3=False, model_id=PDF_CONVERT_MODEL)
                            
                            if markdown_content:
                                # Create markdown filename and save locally
                                base_name = os.path.splitext(filepath)[0]
                                md_filepath = f"{base_name}.md"
                                
                                try:
                                    with open(md_filepath, 'w', encoding='utf-8') as md_file:
                                        md_file.write(markdown_content)
                                    logger.info(f"Saved markdown locally: {md_filepath}")
                                    
                                    # Store both PDF and markdown info for Bedrock processing
                                    new_pdf_files.append({
                                        'path': filepath,
                                        'md_path': md_filepath,
                                        'filename': file.filename,
                                        'timestamped_filename': timestamped_filename,
                                        'has_markdown': True
                                    })
                                    
                                    # Also add markdown content to text context for immediate use
                                    new_text_content += f"\n\n--- Content from {file.filename} (converted from PDF) ---\n{markdown_content}"
                                    
                                except Exception as e:
                                    logger.error(f"Error saving markdown locally: {e}")
                                    # Fall back to PDF-only storage
                                    new_pdf_files.append({
                                        'path': filepath,
                                        'filename': file.filename,
                                        'timestamped_filename': timestamped_filename,
                                        'has_markdown': False
                                    })
                            else:
                                logger.warning(f"Failed to convert PDF to markdown: {file.filename}")
                                # Store PDF info for direct Bedrock processing
                                new_pdf_files.append({
                                    'path': filepath,
                                    'filename': file.filename,
                                    'timestamped_filename': timestamped_filename,
                                    'has_markdown': False
                                })
                            
                            uploaded_files.append({
                                'filename': file.filename,
                                'size': os.path.getsize(filepath),
                                'status': 'success',
                                'type': 'pdf'
                            })
                        else:
                            # Regular text content
                            new_text_content += f"\n\n--- Content from {file.filename} ---\n{content}"
                            uploaded_files.append({
                                'filename': file.filename,
                                'size': os.path.getsize(filepath),
                                'status': 'success',
                                'type': 'text'
                            })
                    else:
                        uploaded_files.append({
                            'filename': file.filename,
                            'status': 'error',
                            'message': 'Could not read file content'
                        })
        
        logger.info(f"Finished processing files. new_pdf_files count: {len(new_pdf_files)}, uploaded_files count: {len(uploaded_files)}")
        logger.info(f"new_pdf_files: {[f['filename'] for f in new_pdf_files]}")
        logger.info(f"uploaded_files: {[f['filename'] for f in uploaded_files]}")
        
        # Process PDFs in parallel if any were queued
        if pdf_conversion_tasks and PDF_CONVERSION_ENABLED:
            logger.info(f"Starting parallel PDF conversion for {len(pdf_conversion_tasks)} files")
            conversion_results = convert_pdfs_parallel(pdf_conversion_tasks)
            
            # Update PDF files with conversion results
            for file_key, result in conversion_results.items():
                if file_key in pdf_files_info:
                    pdf_info = pdf_files_info[file_key]['pdf_info']
                    session_location = pdf_files_info[file_key]['session_location']
                    timestamped_filename = pdf_files_info[file_key]['timestamped_filename']
                    
                    if result['markdown_content']:
                        # Create markdown filename
                        base_name = os.path.splitext(timestamped_filename)[0]
                        md_filename = f"{base_name}.md"
                        md_s3_key = f"{session_location}{md_filename}"
                        
                        # Upload markdown to S3
                        try:
                            s3_client.put_object(
                                Bucket=s3_bucket_name,
                                Key=md_s3_key,
                                Body=result['markdown_content'].encode('utf-8'),
                                ContentType='text/markdown'
                            )
                            logger.info(f"Uploaded markdown to S3: s3://{s3_bucket_name}/{md_s3_key}")
                            
                            # Update PDF info with markdown details
                            pdf_info['md_s3_key'] = md_s3_key
                            pdf_info['has_markdown'] = True
                            
                            # Also add markdown content to text context for immediate use
                            new_text_content += f"\n\n--- Content from {result['filename']} (converted from PDF) ---\n{result['markdown_content']}"
                            
                        except Exception as e:
                            logger.error(f"Error uploading markdown to S3 for {result['filename']}: {e}")
                    else:
                        logger.warning(f"Failed to convert PDF to markdown: {result['filename']} - {result.get('error', 'Unknown error')}")
        else:
            logger.info("No PDF conversion tasks to process")
        
        # Initialize session lists if they don't exist
        if 'uploaded_files' not in session:
            session['uploaded_files'] = []
        if 'pdf_files' not in session:
            session['pdf_files'] = []
        
        logger.info(f"Session state before processing: {len(session['uploaded_files'])} uploaded files, {len(session['pdf_files'])} PDF files")
        logger.info(f"Existing files in session: {[f['filename'] for f in session['uploaded_files']]}")
        
        # Make copies of session lists to avoid concurrent modification issues
        current_uploaded_files = list(session['uploaded_files'])
        current_pdf_files = list(session['pdf_files'])
        
        # Add new PDF files to existing list (avoid duplicates)
        for new_pdf in new_pdf_files:
            existing = next((f for f in current_pdf_files if f['filename'] == new_pdf['filename']), None)
            if not existing:
                current_pdf_files.append(new_pdf)
                logger.info(f"Added new PDF to session: {new_pdf['filename']}")
        
        # Add new files to uploaded files list (avoid duplicates)
        logger.info(f"Before adding new files - session has {len(current_uploaded_files)} files")
        for file_info in uploaded_files:
            if file_info['status'] == 'success':
                existing = next((f for f in current_uploaded_files if f['filename'] == file_info['filename']), None)
                if not existing:
                    current_uploaded_files.append({
                        'filename': file_info['filename'],
                        'size': file_info['size'],
                        'upload_time': datetime.now().isoformat(),
                        'type': file_info.get('type', 'text')
                    })
                    logger.info(f"Added new file to session: {file_info['filename']}")
                else:
                    logger.info(f"File already exists in session: {file_info['filename']}")
        
        # Update session with the modified lists - do this atomically
        session['uploaded_files'] = current_uploaded_files
        session['pdf_files'] = current_pdf_files
        
        logger.info(f"After adding new files - session has {len(current_uploaded_files)} files")
        logger.info(f"Session files: {[f['filename'] for f in current_uploaded_files]}")
        
        # Rebuild complete document context from all text files in session
        total_content = ""
        for file_info in current_uploaded_files:
            if file_info.get('type') != 'pdf':
                if USE_S3_BUCKET:
                    # S3 mode - find file by S3 key pattern
                    try:
                        paginator = s3_client.get_paginator('list_objects_v2')
                        for page in paginator.paginate(Bucket=s3_bucket_name, Prefix=session.get('s3_prefix', '')):
                            if 'Contents' in page:
                                for obj in page['Contents']:
                                    if obj['Key'].endswith(f"_{secure_filename(file_info['filename'])}"):
                                        content = read_file_content(obj['Key'], is_s3=True)
                                        if content and not isinstance(content, dict):
                                            total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                                        break
                    except Exception as e:
                        logger.error(f"Error reading S3 file for context rebuild: {e}")
                else:
                    # Local mode - find file on disk
                    session_folder = session.get('session_folder', '')
                    if os.path.exists(session_folder):
                        for uploaded_file in os.listdir(session_folder):
                            if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                                filepath = os.path.join(session_folder, uploaded_file)
                                content = read_file_content(filepath, is_s3=False)
                                if content and not isinstance(content, dict):
                                    total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                                break
        
        session['document_context'] = f"Context from uploaded documents:{total_content}" if total_content else ""
        session.modified = True
        
        logger.info(f"Session updated - final file count: {len(current_uploaded_files)} uploaded, {len(current_pdf_files)} PDFs")
        
        # Removed redundant log - already logged above
        
        # Calculate total context length including PDFs
        total_context_length = len(total_content)
        pdf_context_note = ""
        if current_pdf_files:
            pdf_count = len(current_pdf_files)
            pdf_context_note = f" + {pdf_count} PDF file(s) processed by Bedrock"
            # Estimate PDF content size for display (rough estimate)
            for pdf_info in current_pdf_files:
                if not USE_S3_BUCKET and 'path' in pdf_info and os.path.exists(pdf_info['path']):
                    pdf_size = os.path.getsize(pdf_info['path'])
                    total_context_length += pdf_size // 4  # Rough estimate: 4 bytes per character
                elif USE_S3_BUCKET and 's3_key' in pdf_info:
                    # For S3, we could get object size but it's optional for display
                    total_context_length += 50000  # Rough estimate for display
        
        context_display = f"{len(total_content)} text chars{pdf_context_note}" if pdf_context_note else f"{total_context_length} chars"
        
        logger.info(f"Returning to frontend - uploaded_files count: {len(current_uploaded_files)}")
        logger.info(f"Files being returned: {[f['filename'] for f in current_uploaded_files]}")
        
        result = {
            'message': f'Successfully processed {len(uploaded_files)} files',
            'files': uploaded_files,
            'context_length': total_context_length,
            'context_display': context_display,
            'uploaded_files': current_uploaded_files,
            'pdf_count': len(current_pdf_files)
        }
        request_end_time = time.time()
        request_duration = request_end_time - request_start_time
        logger.info(f"=== UPLOAD REQUEST END - Duration: {request_duration:.2f}s - Returning result ===")
        return jsonify(result)
        
    except Exception as e:
        logger.error(f"Upload error: {e}")
        return jsonify({'error': 'Upload failed'}), 500

@app.route('/clear')
def clear_session():
    global session_upload_folders
    
    if USE_S3_BUCKET:
        # S3 mode - delete all objects with session prefix
        if 's3_prefix' in session:
            try:
                paginator = s3_client.get_paginator('list_objects_v2')
                for page in paginator.paginate(Bucket=s3_bucket_name, Prefix=session['s3_prefix']):
                    if 'Contents' in page:
                        objects = [{'Key': obj['Key']} for obj in page['Contents']]
                        s3_client.delete_objects(
                            Bucket=s3_bucket_name,
                            Delete={'Objects': objects}
                        )
                logger.info(f"Cleaned up S3 session files with prefix: {session['s3_prefix']}")
            except Exception as e:
                logger.error(f"Error cleaning up S3 session files: {e}")
    else:
        # Local mode - clean up session-specific upload folder
        if 'session_folder' in session:
            session_folder = session['session_folder']
            if os.path.exists(session_folder):
                shutil.rmtree(session_folder)
                logger.info(f"Cleaned up session folder: {session_folder}")
            # Remove from tracking set
            session_upload_folders.discard(session_folder)
    
    session.clear()
    return redirect(url_for('index'))

@app.route('/history')
def get_history():
    return jsonify(session.get('chat_history', []))

@app.route('/files')
def get_files():
    return jsonify(session.get('uploaded_files', []))

@app.route('/remove_file', methods=['POST'])
def remove_file():
    try:
        data = request.get_json()
        filename = data.get('filename')
        
        if not filename:
            return jsonify({'error': 'No filename provided'}), 400
        
        # Remove from uploaded files list
        if 'uploaded_files' in session:
            session['uploaded_files'] = [f for f in session['uploaded_files'] if f['filename'] != filename]
        
        # Remove from PDF files list if it's a PDF and also get markdown info for cleanup
        pdf_to_remove = None
        if 'pdf_files' in session:
            pdf_to_remove = next((f for f in session['pdf_files'] if f['filename'] == filename), None)
            session['pdf_files'] = [f for f in session['pdf_files'] if f['filename'] != filename]
        
        # Remove the actual file from storage (including markdown files)
        if USE_S3_BUCKET:
            # S3 mode - find and delete file from S3
            if 's3_prefix' not in session:
                return jsonify({'error': 'Session S3 prefix not found'}), 400
            
            try:
                # Find and delete both PDF and markdown files in S3
                paginator = s3_client.get_paginator('list_objects_v2')
                for page in paginator.paginate(Bucket=s3_bucket_name, Prefix=session['s3_prefix']):
                    if 'Contents' in page:
                        for obj in page['Contents']:
                            obj_key = obj['Key']
                            secured_filename = secure_filename(filename)
                            
                            # Delete PDF file
                            if obj_key.endswith(f"_{secured_filename}"):
                                s3_client.delete_object(Bucket=s3_bucket_name, Key=obj_key)
                                logger.info(f"Removed file from S3: s3://{s3_bucket_name}/{obj_key}")
                            
                            # Delete corresponding markdown file if it exists
                            elif filename.lower().endswith('.pdf'):
                                base_filename = os.path.splitext(secured_filename)[0]
                                if obj_key.endswith(f"_{base_filename}.md"):
                                    s3_client.delete_object(Bucket=s3_bucket_name, Key=obj_key)
                                    logger.info(f"Removed markdown file from S3: s3://{s3_bucket_name}/{obj_key}")
                                    
                # Also delete specific markdown files if we have the S3 key
                if pdf_to_remove and pdf_to_remove.get('md_s3_key'):
                    try:
                        s3_client.delete_object(Bucket=s3_bucket_name, Key=pdf_to_remove['md_s3_key'])
                        logger.info(f"Removed specific markdown file from S3: s3://{s3_bucket_name}/{pdf_to_remove['md_s3_key']}")
                    except Exception as e:
                        logger.warning(f"Could not remove specific markdown file from S3: {e}")
                        
            except Exception as e:
                logger.error(f"Error removing file from S3: {e}")
        else:
            # Local mode - remove file from disk (including markdown)
            if 'session_folder' not in session or not os.path.exists(session['session_folder']):
                return jsonify({'error': 'Session folder not found'}), 400
            
            session_folder = session['session_folder']
            secured_filename = secure_filename(filename)
            
            for uploaded_file in os.listdir(session_folder):
                file_path = os.path.join(session_folder, uploaded_file)
                
                # Remove PDF file
                if uploaded_file.endswith(f"_{secured_filename}"):
                    try:
                        os.remove(file_path)
                        logger.info(f"Removed file from disk: {file_path}")
                    except Exception as e:
                        logger.error(f"Error removing file {file_path}: {e}")
                
                # Remove corresponding markdown file if it exists
                elif filename.lower().endswith('.pdf'):
                    base_filename = os.path.splitext(secured_filename)[0]
                    if uploaded_file.endswith(f"_{base_filename}.md"):
                        try:
                            os.remove(file_path)
                            logger.info(f"Removed markdown file from disk: {file_path}")
                        except Exception as e:
                            logger.error(f"Error removing markdown file {file_path}: {e}")
            
            # Also remove specific markdown file if we have the path
            if pdf_to_remove and pdf_to_remove.get('md_path') and os.path.exists(pdf_to_remove['md_path']):
                try:
                    os.remove(pdf_to_remove['md_path'])
                    logger.info(f"Removed specific markdown file from disk: {pdf_to_remove['md_path']}")
                except Exception as e:
                    logger.warning(f"Could not remove specific markdown file from disk: {e}")
        
        # Rebuild document context without the removed file
        remaining_files = session.get('uploaded_files', [])
        
        if remaining_files:
            # Re-read remaining text files to rebuild context
            total_content = ""
            for file_info in remaining_files:
                if file_info.get('type') != 'pdf':  # Only process text files
                    if USE_S3_BUCKET:
                        # S3 mode - find file by S3 key pattern
                        try:
                            paginator = s3_client.get_paginator('list_objects_v2')
                            for page in paginator.paginate(Bucket=s3_bucket_name, Prefix=session.get('s3_prefix', '')):
                                if 'Contents' in page:
                                    for obj in page['Contents']:
                                        if obj['Key'].endswith(f"_{secure_filename(file_info['filename'])}"):
                                            content = read_file_content(obj['Key'], is_s3=True)
                                            if content and not isinstance(content, dict):
                                                total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                                            break
                        except Exception as e:
                            logger.error(f"Error reading S3 file for context rebuild: {e}")
                    else:
                        # Local mode - find file on disk
                        session_folder = session.get('session_folder', '')
                        if os.path.exists(session_folder):
                            for uploaded_file in os.listdir(session_folder):
                                if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                                    filepath = os.path.join(session_folder, uploaded_file)
                                    content = read_file_content(filepath, is_s3=False)
                                    if content and not isinstance(content, dict):
                                        total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                                    break
            
            session['document_context'] = f"Context from uploaded documents:{total_content}" if total_content else ""
        else:
            session['document_context'] = ""
        
        session.modified = True
        
        return jsonify({
            'message': f'File {filename} removed successfully',
            'uploaded_files': session.get('uploaded_files', [])
        })
        
    except Exception as e:
        logger.error(f"Remove file error: {e}")
        return jsonify({'error': 'Failed to remove file'}), 500

if __name__ == '__main__':
    # Register signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Register cleanup function for normal exit
    atexit.register(cleanup_resources)

    # Initialize AWS clients
    initialize_aws_clients()

    # Handle S3 bucket creation/loading based on Flask startup mode
    if USE_S3_BUCKET and not os.environ.get('WERKZEUG_RUN_MAIN'):
        # Initial startup - create new bucket
        bucket_created = create_s3_bucket()
        if not bucket_created:
            logger.info("Switched to local filesystem mode (--no-bucket equivalent)")
    elif USE_S3_BUCKET and os.environ.get('WERKZEUG_RUN_MAIN'):
        # Flask debug restart - load existing bucket
        handle_flask_restart_bucket()    

    try:
        logger.info("Starting BedBot application...")
        app.run(debug=True, host='0.0.0.0', port=5000)
    except KeyboardInterrupt:
        # This shouldn't be reached due to signal handler, but kept as fallback
        logger.info("Application interrupted")
    except Exception as e:
        logger.error(f"Application error: {e}")
    finally:
        # Signal handler should handle cleanup, but ensure it runs
        if not cleanup_performed:
            cleanup_resources()

