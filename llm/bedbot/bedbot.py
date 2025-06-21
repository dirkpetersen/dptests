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
import pymupdf4llm
import re
from dotenv import load_dotenv
import signal
import sys
import threading


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
# PDF merging configuration removed - no longer supported
# PDF conversion configuration - temporarily disable to test
PDF_CONVERSION_ENABLED = os.getenv('PDF_CONVERSION', '1').strip().lower() in ('1', 'true', 'yes', 'on')
# Local PDF conversion using fitz instead of Bedrock
PDF_LOCAL_CONVERT = os.getenv('PDF_LOCAL_CONVERT', '0').strip().lower() in ('1', 'true', 'yes', 'on')
# Maximum concurrent Bedrock connections for parallel PDF conversion
BEDROCK_MAX_CONCURRENT = int(os.getenv('BEDROCK_MAXCON', '3'))
# PDF conversion model - use different model for PDF conversion if specified
PDF_CONVERT_MODEL = os.getenv('PDF_CONVERT_MODEL', BEDROCK_MODEL)
# Bedrock timeout configuration
BEDROCK_TIMEOUT = int(os.getenv('BEDROCK_TIMEOUT', '900'))
# PDF chunking configuration removed - no longer supported

if DEBUG_MODE:
    logger.info(f"PDF conversion enabled: {PDF_CONVERSION_ENABLED} (PDF_CONVERSION={os.getenv('PDF_CONVERSION', '1')})")
    logger.info(f"PDF local convert enabled: {PDF_LOCAL_CONVERT} (PDF_LOCAL_CONVERT={os.getenv('PDF_LOCAL_CONVERT', '0')})")
    logger.info(f"Max concurrent Bedrock connections: {BEDROCK_MAX_CONCURRENT} (BEDROCK_MAXCON={os.getenv('BEDROCK_MAXCON', '3')})")
    logger.info(f"PDF conversion model: {PDF_CONVERT_MODEL} (PDF_CONVERT_MODEL={os.getenv('PDF_CONVERT_MODEL', 'default')})")
    logger.info(f"Chat model: {BEDROCK_MODEL} (BEDROCK_MODEL={os.getenv('BEDROCK_MODEL', 'default')}, --model={args.model})")
    # Log parsed models for debugging
    parsed_models = [model.strip() for model in BEDROCK_MODEL.split(',') if model.strip()]
    logger.info(f"Parsed models: {parsed_models} (count: {len(parsed_models)})")
    logger.info(f"Bedrock timeout: {BEDROCK_TIMEOUT}s (BEDROCK_TIMEOUT={os.getenv('BEDROCK_TIMEOUT', '900')})")
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
    
    # Debug logging to understand the environment
    werkzeug_run_main = os.environ.get('WERKZEUG_RUN_MAIN')
    logger.info(f"cleanup_resources called: WERKZEUG_RUN_MAIN={werkzeug_run_main}, PID={os.getpid()}")
    
    # In debug mode with reloader:
    # - Parent process (reloader): WERKZEUG_RUN_MAIN is not set
    # - Child process (actual app): WERKZEUG_RUN_MAIN is set to 'true'
    # We want cleanup to happen in the child process only
    if werkzeug_run_main is None:
        # This is the reloader process in debug mode - don't cleanup
        logger.info("Skipping cleanup - this is the Flask reloader process (parent)")
        return
    
    # Prevent duplicate cleanup
    if cleanup_performed:
        logger.info("Skipping cleanup - already performed")
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
            read_timeout=BEDROCK_TIMEOUT,  # Configurable timeout via BEDROCK_TIMEOUT env var
            connect_timeout=60,  # 1 minute for connection
            retries={
                'max_attempts': 3,
                'mode': 'adaptive'
            }
        )
        
        bedrock_client = session_aws.client('bedrock-runtime', region_name=profile_region, config=config)
        s3_client = session_aws.client('s3', region_name=profile_region, config=config) if USE_S3_BUCKET else None
        
        if DEBUG_MODE:
            logger.info(f"Initialized Bedrock client with profile: {session_aws.profile_name or 'default'}")
            logger.info(f"Using region: {profile_region}")
            logger.info(f"Configured clients with {BEDROCK_TIMEOUT}s read timeout for large document processing")
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


# PDF merging functionality removed - no longer supported

# PDF filename merging helper removed - no longer needed

# PDF chunking functionality removed - no longer supported

# PDF chunk conversion functionality removed - no longer supported

# PDF parallel chunking functionality removed - no longer supported

try:
    from io import BytesIO
    from docling.datamodel.base_models import DocumentStream, InputFormat
    from docling.document_converter import DocumentConverter, PdfFormatOption
    from docling.datamodel.pipeline_options import PdfPipelineOptions
    from docling.backend.pypdfium2_backend import PyPdfiumDocumentBackend
    import ballermann
    DOCLING_AVAILABLE = True
except ImportError:
    DOCLING_AVAILABLE = False
    logger.warning("docling not found - will use fitz text extraction")

def pdf2markdown(pdf_input, is_bytes=False):
    """
    Convert PDF to markdown using docling if available, otherwise fitz text extraction
    Args:
        pdf_input: Path to PDF file or PDF bytes
        is_bytes: True if pdf_input is bytes, False if it's a file path
    Returns:
        Markdown string ready for LLM consumption
    """
    try:
        if DOCLING_AVAILABLE:
            try:
                # PyPdfium without EasyOCR - as per docling example
                pipeline_options = PdfPipelineOptions()
                pipeline_options.do_ocr = False
                pipeline_options.do_table_structure = True
                pipeline_options.table_structure_options.do_cell_matching = False
                
                # Disable GPU usage to avoid CUDA errors
                if hasattr(pipeline_options, 'ocr_options'):
                    pipeline_options.ocr_options.use_gpu = False
                if hasattr(pipeline_options, 'table_structure_options') and hasattr(pipeline_options.table_structure_options, 'use_gpu'):
                    pipeline_options.table_structure_options.use_gpu = False
                
                # Initialize converter with PyPdfium backend (no OCR)
                converter = DocumentConverter(
                    format_options={
                        InputFormat.PDF: PdfFormatOption(
                            pipeline_options=pipeline_options, 
                            backend=PyPdfiumDocumentBackend
                        )
                    }
                )
                
                logger.info("Docling configured with PyPdfium backend: OCR=False, Tables=True, Cell matching=False, GPU=False")
                
                # Use docling for better PDF conversion with table focus
                if is_bytes:
                    # pdf_input is already bytes
                    buf = BytesIO(pdf_input)
                    source = DocumentStream(name="document.pdf", stream=buf)
                else:
                    # Convert file path for docling
                    source = pdf_input
                
                result = converter.convert(source)
                
                # Get markdown content from docling result
                if hasattr(result, 'document') and hasattr(result.document, 'export_to_markdown'):
                    markdown_content = result.document.export_to_markdown()
                    if markdown_content and markdown_content.strip():
                        doc_header = "# Document\n\n*Converted from PDF using docling*\n\n"
                        return doc_header + markdown_content.strip()
                
                # If docling result doesn't have expected structure, fall back to text
                logger.warning("Docling conversion successful but couldn't extract markdown, using fitz fallback")
                
            except Exception as docling_error:
                error_msg = str(docling_error)
                if "meta tensor" in error_msg or "torch" in error_msg.lower():
                    logger.warning(f"Docling PyTorch/GPU error, falling back to fitz: {error_msg}")
                else:
                    logger.warning(f"Docling conversion failed, falling back to fitz: {error_msg}")
                # Continue to fitz fallback
        
        # Fallback to fitz text extraction
        if is_bytes:
            doc = fitz.open(stream=pdf_input, filetype="pdf")
        else:
            doc = fitz.open(pdf_input)
        
        markdown_parts = []
        
        # Process each page with fitz
        for page_num in range(len(doc)):
            page = doc[page_num]
            page_markdown = f"\n## Page {page_num + 1}\n\n"
            
            try:
                # Use basic text extraction
                text = page.get_text("text")
                if text and text.strip():
                    page_markdown += f"{text.strip()}\n\n"
                    
            except Exception as e:
                logger.warning(f"Error processing page {page_num + 1}: {e}")
                page_markdown += "*Error extracting page content*\n\n"
            
            if page_markdown.strip() != f"## Page {page_num + 1}":
                markdown_parts.append(page_markdown)
        
        doc.close()
        
        # Assemble final document
        if not markdown_parts:
            return "# Document\n\n*No content could be extracted from this PDF*"
        
        full_content = "".join(markdown_parts)
        conversion_method = "docling (fallback to fitz)" if DOCLING_AVAILABLE else "fitz"
        doc_header = f"# Document\n\n*Converted from PDF using {conversion_method}*\n\n"
        doc_footer = f"\n\n---\n*Document processed: {len(markdown_parts)} pages*"
        
        return doc_header + full_content + doc_footer
        
    except Exception as e:
        logger.error(f"Error in pdf2markdown conversion: {e}")
        return f"# Document Conversion Error\n\n*Failed to convert PDF: {str(e)}*"

def _convert_table_to_markdown(table_data):
    """Convert table data to markdown format"""
    if not table_data or len(table_data) < 1:
        return None
    
    try:
        markdown_lines = []
        
        # Header row
        header = table_data[0]
        if not header:
            return None
        
        # Clean header cells
        clean_header = [str(cell).strip().replace('\n', ' ') if cell else '' for cell in header]
        markdown_lines.append('| ' + ' | '.join(clean_header) + ' |')
        
        # Separator row
        separator = ['---'] * len(clean_header)
        markdown_lines.append('| ' + ' | '.join(separator) + ' |')
        
        # Data rows
        for row in table_data[1:]:
            if row:
                # Ensure row has same number of columns as header
                clean_row = []
                for i in range(len(clean_header)):
                    if i < len(row) and row[i]:
                        cell_content = str(row[i]).strip().replace('\n', ' ')
                    else:
                        cell_content = ''
                    clean_row.append(cell_content)
                
                markdown_lines.append('| ' + ' | '.join(clean_row) + ' |')
        
        return '\n'.join(markdown_lines)
        
    except Exception as e:
        logger.warning(f"Error converting table to markdown: {e}")
        return None

def convert_pdf_to_markdown_local(pdf_path_or_s3key, is_s3=False, s3_bucket=None):
    """Convert a PDF file to markdown using advanced PDF2Markdown class"""
    doc_name = os.path.basename(pdf_path_or_s3key).replace('.pdf', '')
    
    logger.info(f"Using advanced PDF2Markdown conversion for {doc_name}")
    
    try:
        if is_s3:
            # S3 mode - download PDF first
            if not s3_bucket:
                logger.error("S3 bucket name required for S3 mode")
                return None
            
            try:
                response = s3_client.get_object(Bucket=s3_bucket, Key=pdf_path_or_s3key)
                pdf_bytes = response['Body'].read()
                logger.info(f"Downloaded PDF from S3 for conversion: {doc_name} ({len(pdf_bytes)} bytes)")
                
                # Use the new pdf2markdown function with bytes
                markdown_content = pdf2markdown(pdf_bytes, is_bytes=True)
                logger.info(f"Advanced PDF conversion successful ({len(markdown_content)} characters)")
                return markdown_content
                
            except Exception as download_error:
                logger.error(f"Failed to download PDF from S3: {download_error}")
                return f"# {doc_name}\n\n*PDF conversion failed and could not download from S3*"
        else:
            # Local mode
            if not os.path.exists(pdf_path_or_s3key):
                logger.error(f"PDF file not found: {pdf_path_or_s3key}")
                return f"# {doc_name}\n\n*PDF file not found*"
            
            try:
                # Use the new pdf2markdown function with file path
                markdown_content = pdf2markdown(pdf_path_or_s3key, is_bytes=False)
                logger.info(f"Advanced PDF conversion successful ({len(markdown_content)} characters)")
                return markdown_content
                
            except Exception as conversion_error:
                logger.error(f"Failed to convert PDF: {conversion_error}")
                return f"# {doc_name}\n\n*PDF conversion failed: {str(conversion_error)}*"
            
    except Exception as e:
        logger.error(f"PDF conversion failed: {e}")
        return f"# {doc_name}\n\n*PDF conversion failed completely*"

def convert_pdf_to_markdown_bedrock(pdf_path_or_s3key, is_s3=False, s3_bucket=None, model_id=None):
    """Convert a PDF file to markdown using Bedrock model"""
    if PDF_LOCAL_CONVERT:
        # If local conversion is enabled, use fitz instead
        return convert_pdf_to_markdown_local(pdf_path_or_s3key, is_s3=is_s3, s3_bucket=s3_bucket)
    
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
            conversion_prompt = """Please convert this PDF document to markdown format in full. 

IMPORTANT INSTRUCTIONS:
- Extract all text content from the document
- Preserve the document structure using appropriate markdown headers (# ## ###)
- Convert tables to markdown table format if present
- Include all relevant text content
- Format lists using markdown list syntax
- Preserve paragraph breaks and spacing
- If there are multiple sections, use appropriate header levels
- Do not add any commentary or explanation - just provide the markdown conversion

VERY IMPORTANT:
- Convert the entire document in full even if it is large, do not truncate or skip any sections
- If the document is too large to convert, return a single error message "document too large"
- If the document is empty or contains no text, return a single error message "empty document"
- If the document is not a valid format, return a single error message "invalid format"
"""
            
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
                        logger.info(f"Model {model_id} doesn't support S3Location, falling back to chunked bytes method")
                        # Fall back to downloading and chunking for large PDFs
                        try:
                            response = s3_client.get_object(Bucket=s3_bucket, Key=pdf_path_or_s3key)
                            pdf_bytes = response['Body'].read()
                            
                            logger.info(f"Downloaded PDF from S3: {doc_name} ({len(pdf_bytes)} bytes)")
                            
                            # Process as single document (chunking removed)
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
                
                # Process as single document (chunking removed)
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
    
    # Determine max workers based on conversion method
    max_workers = BEDROCK_MAX_CONCURRENT if not PDF_LOCAL_CONVERT else min(len(pdf_conversion_tasks), 8)
    logger.info(f"Starting parallel PDF conversion for {len(pdf_conversion_tasks)} files (max concurrent: {max_workers}, method: {'local' if PDF_LOCAL_CONVERT else 'bedrock'})")
    
    results = {}
    
    def convert_single_pdf(task):
        """Convert a single PDF - wrapper for thread pool"""
        file_key, s3_key, s3_bucket, filename = task
        try:
            logger.info(f"[Parallel] Starting conversion: {filename}")
            if PDF_LOCAL_CONVERT:
                markdown_content = convert_pdf_to_markdown_local(s3_key, is_s3=True, s3_bucket=s3_bucket)
            else:
                markdown_content = convert_pdf_to_markdown_bedrock(s3_key, is_s3=True, s3_bucket=s3_bucket, model_id=PDF_CONVERT_MODEL)
            logger.info(f"[Parallel] Completed conversion: {filename} ({len(markdown_content) if markdown_content else 0} chars)")
            return file_key, markdown_content, None
        except Exception as e:
            logger.error(f"[Parallel] Failed conversion: {filename} - {e}")
            return file_key, None, str(e)
    
    # Use ThreadPoolExecutor for parallel processing
    try:
        # In debug mode, reduce max_workers to prevent reloader conflicts
        if DEBUG_MODE and max_workers > 2:
            max_workers = 2
            logger.info(f"Reduced max_workers to {max_workers} for debug mode")
            
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_task = {
                executor.submit(convert_single_pdf, task): task 
                for task in pdf_conversion_tasks
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    file_key, markdown_content, error = future.result()
                    
                    results[file_key] = {
                        'markdown_content': markdown_content,
                        'error': error,
                        'filename': task[3]  # filename is the 4th element
                    }
                except Exception as e:
                    logger.error(f"Error getting result for task {task}: {e}")
                    # Add error result for failed task
                    file_key = f"{task[3]}_{len(results)}"  # Generate unique key
                    results[file_key] = {
                        'markdown_content': None,
                        'error': str(e),
                        'filename': task[3]
                    }
    except Exception as e:
        logger.error(f"Critical error in parallel PDF processing: {e}")
        import traceback
        logger.error(f"Full traceback: {traceback.format_exc()}")
        # Return empty results to prevent app crash
        return {}
    
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

def build_pdf_content(pdf_files):
    """Build PDF content for Bedrock messages, preferring markdown when available"""
    content = []
    
    for pdf_info in pdf_files:
        try:
            # Check if we have markdown version available
            has_markdown = pdf_info.get('has_markdown', False)
            doc_name = pdf_info['filename'].replace('.pdf', '')
            doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
            doc_name = re.sub(r'\s+', ' ', doc_name).strip()
            doc_name = doc_name[:50]
            
            if DEBUG_MODE:
                logger.info(f"build_pdf_content: Processing {pdf_info['filename']}, has_markdown={has_markdown}")
                if has_markdown:
                    logger.info(f"build_pdf_content: Markdown keys in pdf_info: {[k for k in pdf_info.keys() if 'md' in k.lower()]}")
            
            if has_markdown:
                # Use markdown content instead of PDF for faster processing
                if USE_S3_BUCKET and 'md_s3_key' in pdf_info:
                    # S3 mode - read markdown from S3
                    try:
                        response = s3_client.get_object(Bucket=s3_bucket_name, Key=pdf_info['md_s3_key'])
                        markdown_content = response['Body'].read().decode('utf-8', errors='ignore')
                        logger.info(f"Using markdown content from S3: {doc_name} ({len(markdown_content)} chars)")
                        
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
                        logger.info(f"Using markdown content from file: {doc_name} ({len(markdown_content)} chars)")
                        
                        content.append({
                            "text": f"\n\n--- Content from {pdf_info['filename']} (PDF converted to markdown) ---\n{markdown_content}"
                        })
                        continue
                    except Exception as e:
                        logger.warning(f"Failed to read markdown from file, falling back to PDF: {e}")
                        has_markdown = False
                else:
                    logger.warning(f"has_markdown=True but no markdown source found for {pdf_info['filename']}")
                    has_markdown = False
            
            # Fall back to PDF processing if no markdown or markdown read failed
            if USE_S3_BUCKET and 's3_key' in pdf_info:
                # S3 mode - try S3 location first, fall back to bytes if model doesn't support it
                s3_key = pdf_info['s3_key']
                
                if 'nooooooooova' in current_model.lower():
                    # Nova models support S3Location
                    try:
                        account_id = session_aws.client('sts').get_caller_identity()['Account']
                        logger.info(f"Adding PDF document from S3 location: {doc_name}")
                        
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
                        response = s3_client.get_object(Bucket=s3_bucket_name, Key=s3_key)
                        pdf_bytes = response['Body'].read()
                        logger.info(f"Adding PDF document (downloaded from S3): {doc_name}")
                        
                        content.append({
                            "document": {
                                "format": "pdf",
                                "name": doc_name,
                                "source": {"bytes": pdf_bytes}
                            }
                        })
                else:
                    # Non-Nova models - use bytes method directly
                    response = s3_client.get_object(Bucket=s3_bucket_name, Key=s3_key)
                    pdf_bytes = response['Body'].read()
                    logger.info(f"Adding PDF document (downloaded from S3): {doc_name}")
                    
                    content.append({
                        "document": {
                            "format": "pdf",
                            "name": doc_name,
                            "source": {"bytes": pdf_bytes}
                        }
                    })
            else:
                # Local mode - read file bytes
                pdf_path = pdf_info['path']
                if not os.path.exists(pdf_path):
                    logger.error(f"PDF file not found: {pdf_path}")
                    continue
                    
                with open(pdf_path, 'rb') as f:
                    pdf_bytes = f.read()
                
                logger.info(f"Adding PDF document: {doc_name}")
                
                content.append({
                    "document": {
                        "format": "pdf",
                        "name": doc_name,
                        "source": {"bytes": pdf_bytes}
                    }
                })
        except Exception as e:
            logger.error(f"Error processing PDF {pdf_info}: {e}")
    
    return content

def send_bedrock_message(messages):
    """Send a single message to Bedrock and return the response"""
    # Use session-selected model or default
    current_model = session.get('selected_model', BEDROCK_MODEL.split(',')[0].strip())
    inference_config = {"maxTokens": 1000, "temperature": 0.5, "topP": 0.9}
    
    logger.info(f"send_bedrock_message using model: {current_model}")
    
    try:
        response = bedrock_client.converse(
            modelId=current_model,
            messages=messages,
            inferenceConfig=inference_config
        )
        return response['output']['message']['content'][0]['text']
    except ClientError as api_error:
        logger.error(f"send_bedrock_message failed with model {current_model}: {api_error}")
        # Handle S3Location fallback
        if 'ValidationException' in str(api_error) and 'S3Location' in str(api_error):
            logger.warning("Model doesn't support S3Location, retrying with bytes method")
            # This would need bytes fallback logic but for now just re-raise
            raise api_error
        else:
            raise api_error

def call_bedrock(prompt, context="", pdf_files=None, conversation_history=None, send_pdfs=False):
    """Call AWS Bedrock model using Converse API with document support and conversation history"""
    if not bedrock_client:
        return "Error: Bedrock client not initialized. Please check your AWS credentials."
    
    # Use session-selected model or default to BEDROCK_MODEL
    current_model = session.get('selected_model', BEDROCK_MODEL.split(',')[0].strip())
    
    # Validate the model is not empty and is in the allowed list
    if not current_model or not current_model.strip():
        current_model = BEDROCK_MODEL.split(',')[0].strip()
        logger.warning(f"Empty model found, using default: {current_model}")
    
    available_models = [model.strip() for model in BEDROCK_MODEL.split(',') if model.strip()]
    if current_model not in available_models:
        current_model = available_models[0]
        logger.warning(f"Invalid model found, using default: {current_model}")
    
    logger.info(f"Using model: {current_model} (from session: {session.get('selected_model', 'not set')})")
    
    # Test if the model is actually available
    try:
        test_response = bedrock_client.converse(
            modelId=current_model,
            messages=[{"role": "user", "content": [{"text": "test"}]}],
            inferenceConfig={"maxTokens": 10, "temperature": 0.1}
        )
        logger.info(f"Model {current_model} is accessible")
    except ClientError as test_error:
        logger.error(f"Model {current_model} test failed: {test_error}")
        if "ValidationException" in str(test_error) and "model identifier is invalid" in str(test_error):
            return f"❌ **Model Error:** The model `{current_model}` is not available in your AWS region or account. Please check:\n\n• Model is available in region `{profile_region}`\n• You have access to this model in AWS Bedrock\n• The model identifier is correct\n\nAvailable models may vary by region and account permissions."
    except Exception as test_error:
        logger.error(f"Model {current_model} test failed with unexpected error: {test_error}")
    
    try:
        messages = []
        
        # Debug logging
        logger.info(f"PDF files provided: {len(pdf_files) if pdf_files else 0}")
        logger.info(f"Send PDFs this call: {send_pdfs}")
        logger.info(f"Conversation history entries: {len(conversation_history) if conversation_history else 0}")
        
        if pdf_files:
            for i, pdf_info in enumerate(pdf_files):
                if USE_S3_BUCKET:
                    logger.info(f"PDF {i+1}: {pdf_info.get('filename', 'unknown')} at s3://{s3_bucket_name}/{pdf_info.get('s3_key', 'unknown')}")
                else:
                    logger.info(f"PDF {i+1}: {pdf_info.get('filename', 'unknown')} at {pdf_info.get('path', 'unknown')}")
        
        # Prepare PDF content if needed
        pdf_content = []
        if send_pdfs and pdf_files:
            logger.info(f"Including PDFs in conversation ({'first time' if not conversation_history else 'model switched'})")
            
            # Build PDF content using existing logic
            pdf_content = build_pdf_content(pdf_files)
                
        # Add conversation history if provided (text-only, no document re-uploads)
        if conversation_history:
            logger.info(f"Adding {len(conversation_history)} conversation history entries as text-only")
            for entry in conversation_history:
                messages.append({"role": "user", "content": [{"text": entry['user']}]})
                messages.append({"role": "assistant", "content": [{"text": entry['bot']}]})
        
        # Create current message prompt  
        if context and not send_pdfs:
            enhanced_prompt = f"{context}\n\n{prompt}"
        else:
            enhanced_prompt = prompt
        
        # Add current user message with PDFs if needed
        current_content = []
        
        # Add PDFs first if this is a re-initialization
        if pdf_content:
            current_content.extend(pdf_content)
            logger.info(f"Added {len(pdf_content)} PDF content items to current message")
        
        # Add the user's text message
        current_content.append({"text": enhanced_prompt})
        
        messages.append({
            "role": "user", 
            "content": current_content
        })
        
        inference_config = {
            "maxTokens": 4000,
            "temperature": 0.7,
            "topP": 0.9
        }
        
        logger.info(f"Sending request to Bedrock with {len(messages)} message entries")
        
        # Debug: Print messages JSON structure if debug mode is enabled
        if DEBUG_MODE and messages:  # Only log if we have messages
            logger.info("=== DEBUG: Converse API Messages Structure ===")
            logger.info(json.dumps(messages, indent=2, default=str))
            logger.info("=== END DEBUG ===")
        
        try:
            response = bedrock_client.converse(
                modelId=current_model,
                messages=messages,
                inferenceConfig=inference_config
            )
            
            return response['output']['message']['content'][0]['text']
            
        except ClientError as api_error:
            # Check if this is an S3Location validation error (only for initial PDF upload)
            if ('ValidationException' in str(api_error) and 'S3Location' in str(api_error) and 
                USE_S3_BUCKET and pdf_files and send_pdfs):
                logger.warning(f"Model {current_model} doesn't support S3Location, retrying with bytes method")
                
                # Rebuild content with bytes method for PDFs (only for initial upload)
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
                    modelId=current_model,
                    messages=retry_messages,
                    inferenceConfig=inference_config
                )
                
                return retry_response['output']['message']['content'][0]['text']
            else:
                # Other API error, re-raise
                raise api_error
        
    except ClientError as e:
        error_message = str(e)
        logger.error(f"Bedrock API error: {e}")
        
        # Check for specific error types that should be shown to the user
        if 'ValidationException' in error_message and 'Input is too long' in error_message:
            return "❌ **Error: Input is too long for the selected model.**\n\nThe uploaded documents and conversation history exceed the model's context limit. Please try:\n\n• Removing some uploaded files\n• Using a model with a larger context window\n• Splitting your request into smaller parts"
        elif 'ValidationException' in error_message:
            return f"❌ **Validation Error:** {error_message}"
        elif 'ThrottlingException' in error_message:
            return "❌ **Rate Limit Exceeded:** Too many requests. Please wait a moment and try again."
        elif 'ServiceUnavailableException' in error_message:
            return "❌ **Service Unavailable:** The Bedrock service is temporarily unavailable. Please try again in a few moments."
        else:
            return f"❌ **Bedrock API Error:** {error_message}"
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return f"❌ **Unexpected Error:** {str(e)}"

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
        session['pdfs_initialized'] = False  # Track if PDFs have been sent to Bedrock
        # Set default model (first in the list)
        default_model = BEDROCK_MODEL.split(',')[0].strip()
        session['selected_model'] = default_model
        logger.info(f"New session created with default model: {default_model}")
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
        chat_history = session.get('chat_history', [])
        pdfs_initialized = session.get('pdfs_initialized', False)
        
        # Debug logging
        logger.info(f"Chat request - PDF files in session: {len(pdf_files)}")
        logger.info(f"Chat request - Chat history length: {len(chat_history)}")
        logger.info(f"Chat request - PDFs initialized: {pdfs_initialized}")
        
        # Determine if we need to send PDFs (first time with PDFs or PDFs not yet initialized)
        send_pdfs = pdf_files and not pdfs_initialized
        
        logger.info(f"PDF decision: pdf_files={len(pdf_files) if pdf_files else 0}, pdfs_initialized={pdfs_initialized}, send_pdfs={send_pdfs}")
        
        # Build conversation history for Bedrock (exclude current session initialization if it exists)
        conversation_history = []
        if chat_history:
            # Filter out any PDF initialization messages (they contain "confirm" language)
            # but only skip the very first message if it's a PDF initialization
            for i, entry in enumerate(chat_history):
                user_msg = entry.get('user', '').lower()
                # Skip only the first message if it's clearly a PDF initialization
                if i == 0 and ('confirm' in user_msg and 'access' in user_msg and ('document' in user_msg or 'uploaded' in user_msg)):
                    logger.info("Skipping PDF initialization message from conversation history")
                    continue
                
                conversation_history.append({
                    'user': entry['user'],
                    'bot': entry['bot']  # Use original markdown, not HTML
                })
        
        logger.info(f"Conversation history entries to send: {len(conversation_history)}")
        logger.info(f"Will send PDFs this request: {send_pdfs}")
        
        # Call Bedrock with conversation history
        bot_response = call_bedrock(
            prompt=user_message,
            context=context,
            pdf_files=pdf_files,
            conversation_history=conversation_history,
            send_pdfs=send_pdfs
        )
        
        # Check if bot_response contains an error (starts with ❌)
        if bot_response.startswith('❌'):
            # Don't mark PDFs as initialized if there was an error
            # Don't store error messages in chat history
            bot_response_html = convert_markdown_to_html(bot_response)
            return jsonify({
                'response': bot_response_html,
                'timestamp': datetime.now().isoformat(),
                'error': True  # Flag to help frontend handle errors differently
            })
        
        # Mark PDFs as initialized if they were sent
        if send_pdfs:
            session['pdfs_initialized'] = True
            logger.info("Marked PDFs as initialized in session")
        
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
        
        # Process all files individually (PDF merging removed)
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
                
                if file_size > MAX_FILE_SIZE:
                    logger.info(f'File "{file.filename}" is {file_size / (1024*1024):.1f} MB')
                
                filename = secure_filename(file.filename)
                # Add timestamp to avoid conflicts
                timestamped_filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{filename}"
                
                if USE_S3_BUCKET:
                    # S3 mode - upload to S3
                    s3_key = f"{session_location}{timestamped_filename}"
                    
                    try:
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
                                        conversion_method = "local" if PDF_LOCAL_CONVERT else "Bedrock"
                                        logger.info(f"Added PDF to parallel conversion queue ({conversion_method}): {file.filename}")
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
                    
                    # Save original file
                    file.save(filepath)
                    
                    # Handle file content based on type
                    content = read_file_content(filepath, is_s3=False)
                    if content:
                        if isinstance(content, dict) and content.get('type') == 'pdf':
                            # Convert PDF to markdown
                            conversion_method = "local" if PDF_LOCAL_CONVERT else "Bedrock"
                            logger.info(f"Converting PDF to markdown ({conversion_method}): {file.filename}")
                            if PDF_LOCAL_CONVERT:
                                markdown_content = convert_pdf_to_markdown_local(filepath, is_s3=False)
                            else:
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
        pdf_added = False
        for new_pdf in new_pdf_files:
            existing = next((f for f in current_pdf_files if f['filename'] == new_pdf['filename']), None)
            if not existing:
                current_pdf_files.append(new_pdf)
                logger.info(f"Added new PDF to session: {new_pdf['filename']}")
                pdf_added = True
        
        # If new PDFs were added, reset the initialization state so they get sent to Bedrock
        if pdf_added:
            session['pdfs_initialized'] = False
            logger.info("Reset PDF initialization state due to new PDF uploads")
        
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
        
        # Calculate total content size sent to LLM (markdown from PDFs + text files)
        total_llm_content_bytes = 0
        
        # Add markdown content from PDFs
        for pdf_info in current_pdf_files:
            if pdf_info.get('has_markdown', False):
                if USE_S3_BUCKET and 'md_s3_key' in pdf_info:
                    # S3 mode - get markdown size from S3
                    try:
                        response = s3_client.head_object(Bucket=s3_bucket_name, Key=pdf_info['md_s3_key'])
                        total_llm_content_bytes += response['ContentLength']
                    except Exception as e:
                        logger.warning(f"Could not get size of markdown file from S3: {e}")
                elif 'md_path' in pdf_info and os.path.exists(pdf_info['md_path']):
                    # Local mode - get markdown size from file
                    try:
                        total_llm_content_bytes += os.path.getsize(pdf_info['md_path'])
                    except Exception as e:
                        logger.warning(f"Could not get size of markdown file: {e}")
        
        # Add text content from regular text files
        for file_info in current_uploaded_files:
            if file_info.get('type') != 'pdf':  # Text files
                if USE_S3_BUCKET:
                    # S3 mode - find file by S3 key pattern and get its size
                    try:
                        paginator = s3_client.get_paginator('list_objects_v2')
                        for page in paginator.paginate(Bucket=s3_bucket_name, Prefix=session.get('s3_prefix', '')):
                            if 'Contents' in page:
                                for obj in page['Contents']:
                                    if obj['Key'].endswith(f"_{secure_filename(file_info['filename'])}"):
                                        total_llm_content_bytes += obj['Size']
                                        break
                    except Exception as e:
                        logger.warning(f"Could not get size of text file from S3: {e}")
                else:
                    # Local mode - find file on disk and get its size
                    session_folder = session.get('session_folder', '')
                    if os.path.exists(session_folder):
                        for uploaded_file in os.listdir(session_folder):
                            if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                                filepath = os.path.join(session_folder, uploaded_file)
                                try:
                                    total_llm_content_bytes += os.path.getsize(filepath)
                                except Exception as e:
                                    logger.warning(f"Could not get size of text file: {e}")
                                break
        
        # Convert to MiB and format with 3 decimal places
        total_llm_content_mib = total_llm_content_bytes / (1024 * 1024)
        llm_content_size_display = f"{total_llm_content_mib:.3f} MiB" if total_llm_content_bytes > 0 else "0.000 MiB"
        
        result = {
            'message': f'Successfully processed {len(uploaded_files)} files',
            'files': uploaded_files,
            'context_length': total_context_length,
            'context_display': context_display,
            'uploaded_files': current_uploaded_files,
            'pdf_count': len(current_pdf_files),
            'markdown_size_display': llm_content_size_display  # Keep same key name for frontend compatibility
        }
        request_end_time = time.time()
        request_duration = request_end_time - request_start_time
        logger.info(f"=== UPLOAD REQUEST END - Duration: {request_duration:.2f}s - Returning result ===")
        return jsonify(result)
        
    except Exception as e:
        import traceback
        logger.error(f"Upload error: {e}")
        logger.error(f"Full traceback: {traceback.format_exc()}")
        return jsonify({'error': f'Upload failed: {str(e)}'}), 500

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

@app.route('/models')
def get_models():
    """Return available Bedrock models from BEDROCK_MODEL configuration"""
    try:
        # Split BEDROCK_MODEL by comma and strip whitespace
        models = [model.strip() for model in BEDROCK_MODEL.split(',') if model.strip()]
        
        # Return models with exact names from config
        model_list = []
        for model in models:
            model_list.append({
                'id': model,
                'name': model  # Use exact model name instead of display name
            })
        
        return jsonify({
            'models': model_list,
            'current': models[0] if models else None  # First model is default
        })
        
    except Exception as e:
        logger.error(f"Error getting models: {e}")
        return jsonify({'models': [], 'current': None})

@app.route('/test_models')
def test_models():
    """Test which models are actually available in AWS Bedrock"""
    try:
        models = [model.strip() for model in BEDROCK_MODEL.split(',') if model.strip()]
        test_results = []
        
        for model in models:
            try:
                # Test each model with a simple request
                test_response = bedrock_client.converse(
                    modelId=model,
                    messages=[{"role": "user", "content": [{"text": "test"}]}],
                    inferenceConfig={"maxTokens": 10, "temperature": 0.1}
                )
                test_results.append({
                    'model': model,
                    'status': 'available',
                    'error': None
                })
                logger.info(f"Model test PASSED: {model}")
            except Exception as e:
                test_results.append({
                    'model': model,
                    'status': 'unavailable',
                    'error': str(e)
                })
                logger.error(f"Model test FAILED: {model} - {e}")
        
        return jsonify({
            'region': profile_region,
            'tests': test_results,
            'summary': {
                'total': len(models),
                'available': len([r for r in test_results if r['status'] == 'available']),
                'unavailable': len([r for r in test_results if r['status'] == 'unavailable'])
            }
        })
        
    except Exception as e:
        logger.error(f"Error testing models: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/change_model', methods=['POST'])
def change_model():
    """Change the active Bedrock model for the session"""
    try:
        data = request.get_json()
        new_model = data.get('model', '').strip()
        
        if not new_model:
            return jsonify({'error': 'No model specified'}), 400
        
        # Validate that the model is in the allowed list
        available_models = [model.strip() for model in BEDROCK_MODEL.split(',') if model.strip()]
        if new_model not in available_models:
            return jsonify({'error': 'Model not available'}), 400
        
        # Store the selected model in the session
        session['selected_model'] = new_model
        
        # Reset only the Bedrock conversation state, keep chat history visible
        session['pdfs_initialized'] = False  # Reset PDF initialization so they get re-sent to new model
        session.modified = True
        
        logger.info(f"Model changed to: {new_model}, PDFs will be re-initialized with new model")
        return jsonify({'success': True, 'model': new_model, 'model_changed': True})
        
    except Exception as e:
        logger.error(f"Error changing model: {e}")
        return jsonify({'error': 'Failed to change model'}), 500

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
            
            # If we removed a PDF and there are still PDFs remaining, reset initialization
            # If no PDFs remain, we can keep the initialization state as is (no PDFs to initialize)
            if pdf_to_remove and session['pdf_files']:
                session['pdfs_initialized'] = False
                logger.info("Reset PDF initialization state due to PDF file removal (PDFs still remain)")
            elif pdf_to_remove and not session['pdf_files']:
                # All PDFs removed - reset state for cleanliness
                session['pdfs_initialized'] = False
                logger.info("Reset PDF initialization state - no PDFs remain")
        
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
    # Only register signal handlers and cleanup in the actual app process
    # In debug mode, this prevents the reloader process from interfering
    if os.environ.get('WERKZEUG_RUN_MAIN') is not None or not args.debug:
        # Register signal handlers for graceful shutdown
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
        
        # Register cleanup function for normal exit
        atexit.register(cleanup_resources)
        logger.info(f"Registered signal handlers and cleanup (debug={args.debug}, WERKZEUG_RUN_MAIN={os.environ.get('WERKZEUG_RUN_MAIN')})")

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

