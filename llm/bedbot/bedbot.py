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

# Create temporary directories
temp_upload_dir = tempfile.mkdtemp(prefix='bedbot_uploads_')
temp_session_dir = tempfile.mkdtemp(prefix='bedbot_sessions_')
app.config['UPLOAD_FOLDER'] = temp_upload_dir
app.config['SESSION_FILE_DIR'] = temp_session_dir

# Initialize Flask-Session
Session(app)

logger.info(f"Created temporary upload directory: {temp_upload_dir}")
logger.info(f"Created temporary session directory: {temp_session_dir}")

# Register cleanup function to remove temp directories on shutdown
def cleanup_temp_dirs():
    if os.path.exists(temp_upload_dir):
        shutil.rmtree(temp_upload_dir)
        logger.info(f"Cleaned up temporary upload directory: {temp_upload_dir}")
    if os.path.exists(temp_session_dir):
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
    """Read content from uploaded file"""
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error reading file {filepath}: {e}")
        return None

def call_bedrock_nova(prompt, context=""):
    """Call AWS Bedrock Nova Lite model using Converse API"""
    if not bedrock_client:
        return "Error: Bedrock client not initialized. Please check your AWS credentials."
    
    try:
        # Combine context and prompt
        full_prompt = f"{context}\n\n{prompt}" if context else prompt
        
        messages = [
            {
                "role": "user",
                "content": [
                    {
                        "text": full_prompt
                    }
                ]
            }
        ]
        
        inference_config = {
            "maxTokens": 4000,
            "temperature": 0.7,
            "topP": 0.9
        }
        
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

@app.route('/')
def index():
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
        session['chat_history'] = []
        session['document_context'] = ""
        session['uploaded_files'] = []
    return render_template('index.html')

@app.route('/chat', methods=['POST'])
def chat():
    try:
        data = request.get_json()
        user_message = data.get('message', '').strip()
        
        if not user_message:
            return jsonify({'error': 'Empty message'}), 400
        
        # Get document context if available
        context = session.get('document_context', '')
        
        # Call Bedrock Nova
        bot_response = call_bedrock_nova(user_message, context)
        
        # Store in chat history
        if 'chat_history' not in session:
            session['chat_history'] = []
        
        session['chat_history'].append({
            'user': user_message,
            'bot': bot_response,
            'timestamp': datetime.now().isoformat()
        })
        
        session.modified = True
        
        return jsonify({
            'response': bot_response,
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
        total_content = ""
        
        for file in files:
            if file and file.filename and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                # Add timestamp to avoid conflicts
                filename = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{filename}"
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                
                file.save(filepath)
                
                # Read file content
                content = read_file_content(filepath)
                if content:
                    total_content += f"\n\n--- Content from {file.filename} ---\n{content}"
                    uploaded_files.append({
                        'filename': file.filename,
                        'size': os.path.getsize(filepath),
                        'status': 'success'
                    })
                else:
                    uploaded_files.append({
                        'filename': file.filename,
                        'status': 'error',
                        'message': 'Could not read file content'
                    })
        
        # Store document context and file list in session
        if total_content:
            session['document_context'] = f"Context from uploaded documents:{total_content}"
            
            # Add to uploaded files list (avoid duplicates)
            if 'uploaded_files' not in session:
                session['uploaded_files'] = []
            
            for file_info in uploaded_files:
                if file_info['status'] == 'success':
                    # Check if file already exists
                    existing = next((f for f in session['uploaded_files'] if f['filename'] == file_info['filename']), None)
                    if not existing:
                        session['uploaded_files'].append({
                            'filename': file_info['filename'],
                            'size': file_info['size'],
                            'upload_time': datetime.now().isoformat()
                        })
            
            session.modified = True
        
        return jsonify({
            'message': f'Successfully processed {len(uploaded_files)} files',
            'files': uploaded_files,
            'context_length': len(total_content),
            'uploaded_files': session.get('uploaded_files', [])
        })
        
    except Exception as e:
        logger.error(f"Upload error: {e}")
        return jsonify({'error': 'Upload failed'}), 500

@app.route('/clear')
def clear_session():
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
        
        # Rebuild document context without the removed file
        remaining_files = session.get('uploaded_files', [])
        if remaining_files:
            # Re-read remaining files to rebuild context
            total_content = ""
            for file_info in remaining_files:
                # Find the actual file on disk (with timestamp prefix)
                for uploaded_file in os.listdir(app.config['UPLOAD_FOLDER']):
                    if uploaded_file.endswith(f"_{secure_filename(file_info['filename'])}"):
                        filepath = os.path.join(app.config['UPLOAD_FOLDER'], uploaded_file)
                        content = read_file_content(filepath)
                        if content:
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
    app.run(debug=True, host='0.0.0.0', port=5000)

