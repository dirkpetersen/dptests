#! /usr/bin/env python3

import os
import json
import uuid
import tempfile
import shutil
import atexit
from datetime import datetime
from flask import Flask, render_template, request, jsonify, session, redirect, url_for
from flask_session import Session
from werkzeug.utils import secure_filename
import boto3
from botocore.exceptions import ClientError
import logging
import markdown

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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

# Track all session upload folders for cleanup
session_upload_folders = set()

# Initialize Flask-Session (will create session dir when first needed)
Session(app)

# Register cleanup function to remove temp directories on shutdown
def cleanup_temp_dirs():
    global temp_session_dir, session_upload_folders
    
    # Clean up all session upload folders
    for folder in list(session_upload_folders):
        if os.path.exists(folder):
            shutil.rmtree(folder)
            logger.info(f"Cleaned up session upload folder: {folder}")
    
    # Clean up Flask session directory
    if temp_session_dir and os.path.exists(temp_session_dir):
        shutil.rmtree(temp_session_dir)
        logger.info(f"Cleaned up temporary session directory: {temp_session_dir}")

atexit.register(cleanup_temp_dirs)

# AWS Bedrock client - uses current AWS profile
try:
    # Create a session that uses the current AWS profile
    session_aws = boto3.Session()
    
    # Get the region from the profile configuration
    profile_region = session_aws.region_name
    if not profile_region:
        # Fallback to default region if not set in profile
        profile_region = 'us-east-1'
        logger.warning("No region found in AWS profile, defaulting to us-east-1")
    
    bedrock_client = session_aws.client('bedrock-runtime', region_name=profile_region)
    logger.info(f"Initialized Bedrock client with profile: {session_aws.profile_name or 'default'}")
    logger.info(f"Using region: {profile_region}")
except Exception as e:
    logger.error(f"Failed to initialize Bedrock client: {e}")
    bedrock_client = None

# Allowed file extensions
ALLOWED_EXTENSIONS = {'txt', 'pdf', 'doc', 'docx', 'md', 'json', 'csv'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def read_file_content(filepath):
    """Read content from uploaded file or return file info for PDF processing"""
    try:
        filename = os.path.basename(filepath).lower()
        
        if filename.endswith('.pdf'):
            # For PDFs, return file info for Nova document processing
            return {
                'type': 'pdf',
                'path': filepath,
                'filename': os.path.basename(filepath)
            }
        else:
            # For text files, read content normally
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                return f.read()
    except Exception as e:
        logger.error(f"Error reading file {filepath}: {e}")
        return None

def call_bedrock_nova(prompt, context="", pdf_files=None):
    """Call AWS Bedrock Nova Lite model using Converse API with document support"""
    if not bedrock_client:
        return "Error: Bedrock client not initialized. Please check your AWS credentials."
    
    try:
        content = []
        
        # Debug logging
        logger.info(f"PDF files provided: {len(pdf_files) if pdf_files else 0}")
        if pdf_files:
            for i, pdf_info in enumerate(pdf_files):
                logger.info(f"PDF {i+1}: {pdf_info.get('filename', 'unknown')} at {pdf_info.get('path', 'unknown')}")
        
        # Add PDF documents if provided
        if pdf_files:
            for pdf_info in pdf_files:
                try:
                    pdf_path = pdf_info['path']
                    if not os.path.exists(pdf_path):
                        logger.error(f"PDF file not found: {pdf_path}")
                        continue
                        
                    with open(pdf_path, 'rb') as f:
                        pdf_bytes = f.read()
                    
                    # Sanitize document name for Nova - only alphanumeric, whitespace, hyphens, parentheses, square brackets
                    doc_name = pdf_info['filename'].replace('.pdf', '')
                    # Replace periods and other invalid characters with underscores
                    import re
                    doc_name = re.sub(r'[^a-zA-Z0-9\s\-\(\)\[\]]', '_', doc_name)
                    # Replace multiple consecutive whitespace with single space
                    doc_name = re.sub(r'\s+', ' ', doc_name).strip()
                    # Limit length
                    doc_name = doc_name[:50]
                    
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
                    logger.error(f"Error reading PDF {pdf_info['path']}: {e}")
        
        # Add text context if available and no PDFs
        if context and not pdf_files:
            enhanced_prompt = f"{context}\n\n{prompt}"
        else:
            # Create enhanced prompt for document analysis
            if pdf_files and len(pdf_files) > 1:
                # When multiple PDFs are provided, give clear instructions for comparison
                enhanced_prompt = f"""I have provided {len(pdf_files)} documents for analysis. Please carefully analyze ALL the provided documents and answer the following question:

{prompt}

IMPORTANT INSTRUCTIONS:
- You have access to {len(pdf_files)} PDF documents - please read and analyze ALL of them
- Compare and analyze the content across all provided documents
- If one document contains requirements or criteria, evaluate the other documents (CVs/resumes) against those criteria
- Provide specific names, examples, and evidence from the documents to support your analysis
- Be thorough and analytical in your comparison
- When asked for the option or to compare options presented in different documents, provide specific and detailed comparisons"""
            elif pdf_files and len(pdf_files) == 1:
                enhanced_prompt = f"""I have provided 1 document for analysis. Please carefully analyze the provided document and answer the following question:

{prompt}

Please provide specific information and examples from the document to support your response."""
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
        
        logger.info(f"Sending request to Nova with {len(content)} content items")
        
        response = bedrock_client.converse(
            modelId='us.amazon.nova-lite-v1:0',
            messages=messages,
            inferenceConfig=inference_config
        )
        
        return response['output']['message']['content'][0]['text']
        
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

def get_session_upload_folder():
    """Get or create session-specific upload folder"""
    global temp_session_dir, session_upload_folders
    
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    
    # Create Flask session directory if it doesn't exist
    if not temp_session_dir:
        temp_session_dir = tempfile.mkdtemp(prefix='bedbot_sessions_')
        os.chmod(temp_session_dir, 0o700)  # Only owner can read/write/execute
        app.config['SESSION_FILE_DIR'] = temp_session_dir
        logger.info(f"Created temporary session directory: {temp_session_dir}")
    
    # Create a unique temporary directory for this session with restricted permissions
    session_folder = tempfile.mkdtemp(prefix=f'bedbot_session_{session["session_id"]}_')
    os.chmod(session_folder, 0o700)  # Only owner can read/write/execute
    
    # Track this folder for cleanup on shutdown
    session_upload_folders.add(session_folder)
    
    # Store the session folder path in the session for cleanup
    session['session_folder'] = session_folder
    session.modified = True
    
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
        # Create session folder
        get_session_upload_folder()
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
        bot_response = call_bedrock_nova(user_message, context, pdf_files)
        
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

@app.route('/upload', methods=['POST'])
def upload_files():
    try:
        if 'files' not in request.files:
            return jsonify({'error': 'No files provided'}), 400
        
        files = request.files.getlist('files')
        uploaded_files = []
        new_text_content = ""
        new_pdf_files = []
        
        # Get or create session-specific upload folder
        if 'session_folder' in session and os.path.exists(session['session_folder']):
            session_folder = session['session_folder']
        else:
            session_folder = get_session_upload_folder()
        
        for file in files:
            if file and file.filename and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                # Add timestamp to avoid conflicts
                timestamped_filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{filename}"
                filepath = os.path.join(session_folder, timestamped_filename)
                
                file.save(filepath)
                
                # Handle file content based on type
                content = read_file_content(filepath)
                if content:
                    if isinstance(content, dict) and content.get('type') == 'pdf':
                        # Store PDF file info for Nova processing
                        new_pdf_files.append({
                            'path': filepath,
                            'filename': file.filename,
                            'timestamped_filename': timestamped_filename
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
        
        # Initialize session lists if they don't exist
        if 'uploaded_files' not in session:
            session['uploaded_files'] = []
        if 'pdf_files' not in session:
            session['pdf_files'] = []
        
        # Add new PDF files to existing list (avoid duplicates)
        for new_pdf in new_pdf_files:
            existing = next((f for f in session['pdf_files'] if f['filename'] == new_pdf['filename']), None)
            if not existing:
                session['pdf_files'].append(new_pdf)
        
        # Add new files to uploaded files list (avoid duplicates)
        for file_info in uploaded_files:
            if file_info['status'] == 'success':
                existing = next((f for f in session['uploaded_files'] if f['filename'] == file_info['filename']), None)
                if not existing:
                    session['uploaded_files'].append({
                        'filename': file_info['filename'],
                        'size': file_info['size'],
                        'upload_time': datetime.now().isoformat(),
                        'type': file_info.get('type', 'text')
                    })
        
        # Rebuild complete document context from all text files in session
        total_content = ""
        for file_info in session['uploaded_files']:
            if file_info.get('type') != 'pdf':
                # Find the actual file on disk
                for uploaded_file in os.listdir(session_folder):
                    if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                        filepath = os.path.join(session_folder, uploaded_file)
                        content = read_file_content(filepath)
                        if content and not isinstance(content, dict):
                            total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                        break
        
        session['document_context'] = f"Context from uploaded documents:{total_content}" if total_content else ""
        session.modified = True
        
        logger.info(f"Session now has {len(session['uploaded_files'])} total files ({len(session['pdf_files'])} PDFs)")
        
        # Calculate total context length including PDFs
        total_context_length = len(total_content)
        pdf_context_note = ""
        if session.get('pdf_files'):
            pdf_count = len(session['pdf_files'])
            pdf_context_note = f" + {pdf_count} PDF file(s) processed by Nova"
            # Estimate PDF content size for display (rough estimate)
            for pdf_info in session['pdf_files']:
                if os.path.exists(pdf_info['path']):
                    pdf_size = os.path.getsize(pdf_info['path'])
                    total_context_length += pdf_size // 4  # Rough estimate: 4 bytes per character
        
        context_display = f"{len(total_content)} text chars{pdf_context_note}" if pdf_context_note else f"{total_context_length} chars"
        
        return jsonify({
            'message': f'Successfully processed {len(uploaded_files)} files',
            'files': uploaded_files,
            'context_length': total_context_length,
            'context_display': context_display,
            'uploaded_files': session.get('uploaded_files', []),
            'pdf_count': len(session.get('pdf_files', []))
        })
        
    except Exception as e:
        logger.error(f"Upload error: {e}")
        return jsonify({'error': 'Upload failed'}), 500

@app.route('/clear')
def clear_session():
    global session_upload_folders
    
    # Clean up session-specific upload folder
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
        
        # Remove from PDF files list if it's a PDF
        if 'pdf_files' in session:
            session['pdf_files'] = [f for f in session['pdf_files'] if f['filename'] != filename]
        
        # Get session folder
        if 'session_folder' not in session or not os.path.exists(session['session_folder']):
            return jsonify({'error': 'Session folder not found'}), 400
        
        session_folder = session['session_folder']
        
        # Rebuild document context without the removed file
        remaining_files = session.get('uploaded_files', [])
        
        if remaining_files:
            # Re-read remaining text files to rebuild context
            total_content = ""
            for file_info in remaining_files:
                if file_info.get('type') != 'pdf':  # Only process text files
                    # Find the actual file on disk (with timestamp prefix)
                    for uploaded_file in os.listdir(session_folder):
                        if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                            filepath = os.path.join(session_folder, uploaded_file)
                            content = read_file_content(filepath)
                            if content and not isinstance(content, dict):
                                total_content += f"\n\n--- Content from {file_info['filename']} ---\n{content}"
                            break
            
            session['document_context'] = f"Context from uploaded documents:{total_content}" if total_content else ""
        else:
            session['document_context'] = ""
        
        # Also remove the actual file from disk
        for uploaded_file in os.listdir(session_folder):
            if uploaded_file.endswith(f"_{secure_filename(filename)}"):
                file_path = os.path.join(session_folder, uploaded_file)
                try:
                    os.remove(file_path)
                    logger.info(f"Removed file from disk: {file_path}")
                except Exception as e:
                    logger.error(f"Error removing file {file_path}: {e}")
                break
        
        session.modified = True
        
        return jsonify({
            'message': f'File {filename} removed successfully',
            'uploaded_files': session.get('uploaded_files', [])
        })
        
    except Exception as e:
        logger.error(f"Remove file error: {e}")
        return jsonify({'error': 'Failed to remove file'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

