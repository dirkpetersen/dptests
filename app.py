import os
import json
import uuid
from datetime import datetime
from flask import Flask, render_template, request, jsonify, session, redirect, url_for
from flask_socketio import SocketIO, emit
from werkzeug.utils import secure_filename
import boto3
from botocore.exceptions import ClientError
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'your-secret-key-here')
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024  # 1GB max file size
app.config['UPLOAD_FOLDER'] = 'uploads'

# Ensure upload directory exists
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

# Initialize SocketIO for real-time communication
socketio = SocketIO(app, cors_allowed_origins="*")

# AWS Bedrock client
try:
    bedrock_client = boto3.client(
        'bedrock-runtime',
        region_name=os.environ.get('AWS_REGION', 'us-east-1')
    )
except Exception as e:
    logger.error(f"Failed to initialize Bedrock client: {e}")
    bedrock_client = None

# Store chat sessions in memory (use Redis/DB in production)
chat_sessions = {}

ALLOWED_EXTENSIONS = {'txt', 'pdf', 'doc', 'docx', 'md', 'json', 'csv'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def get_session_id():
    if 'session_id' not in session:
        session['session_id'] = str(uuid.uuid4())
    return session['session_id']

def call_bedrock_nova(prompt, context=""):
    """Call AWS Bedrock Nova Lite model"""
    if not bedrock_client:
        return "AWS Bedrock client not initialized. Please check your AWS credentials."
    
    try:
        # Prepare the request body for Nova Lite
        body = {
            "inputText": f"{context}\n\nUser: {prompt}\n\nAssistant:",
            "textGenerationConfig": {
                "maxTokenCount": 1000,
                "temperature": 0.7,
                "topP": 0.9
            }
        }
        
        response = bedrock_client.invoke_model(
            modelId="us.amazon.nova-lite-v1:0",
            body=json.dumps(body),
            contentType="application/json",
            accept="application/json"
        )
        
        response_body = json.loads(response['body'].read())
        return response_body.get('results', [{}])[0].get('outputText', 'No response generated.')
        
    except ClientError as e:
        logger.error(f"Bedrock API error: {e}")
        return f"Error calling Bedrock: {str(e)}"
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return f"Unexpected error: {str(e)}"

@app.route('/')
def index():
    session_id = get_session_id()
    if session_id not in chat_sessions:
        chat_sessions[session_id] = {
            'messages': [],
            'documents': [],
            'context': ""
        }
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    session_id = get_session_id()
    
    if 'files' not in request.files:
        return jsonify({'error': 'No files provided'}), 400
    
    files = request.files.getlist('files')
    uploaded_files = []
    
    for file in files:
        if file and file.filename and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            # Add timestamp to avoid conflicts
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_')
            filename = timestamp + filename
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            
            try:
                file.save(filepath)
                
                # Read file content for context
                with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()[:10000]  # Limit to first 10k chars
                
                file_info = {
                    'filename': file.filename,
                    'filepath': filepath,
                    'content': content,
                    'uploaded_at': datetime.now().isoformat()
                }
                
                chat_sessions[session_id]['documents'].append(file_info)
                uploaded_files.append(file.filename)
                
                # Update context with document content
                chat_sessions[session_id]['context'] += f"\n\nDocument: {file.filename}\nContent: {content[:1000]}..."
                
            except Exception as e:
                logger.error(f"Error processing file {file.filename}: {e}")
                return jsonify({'error': f'Error processing {file.filename}'}), 500
    
    return jsonify({
        'message': f'Successfully uploaded {len(uploaded_files)} files',
        'files': uploaded_files
    })

@socketio.on('send_message')
def handle_message(data):
    session_id = get_session_id()
    user_message = data.get('message', '').strip()
    
    if not user_message:
        return
    
    # Initialize session if needed
    if session_id not in chat_sessions:
        chat_sessions[session_id] = {
            'messages': [],
            'documents': [],
            'context': ""
        }
    
    # Add user message to session
    user_msg = {
        'type': 'user',
        'content': user_message,
        'timestamp': datetime.now().isoformat()
    }
    chat_sessions[session_id]['messages'].append(user_msg)
    
    # Emit user message to client
    emit('message_response', user_msg)
    
    # Get AI response
    context = chat_sessions[session_id]['context']
    ai_response = call_bedrock_nova(user_message, context)
    
    # Add AI response to session
    ai_msg = {
        'type': 'assistant',
        'content': ai_response,
        'timestamp': datetime.now().isoformat()
    }
    chat_sessions[session_id]['messages'].append(ai_msg)
    
    # Emit AI response to client
    emit('message_response', ai_msg)

@socketio.on('connect')
def handle_connect():
    session_id = get_session_id()
    logger.info(f'Client connected: {session_id}')
    
    # Send chat history if exists
    if session_id in chat_sessions:
        for message in chat_sessions[session_id]['messages']:
            emit('message_response', message)

@socketio.on('disconnect')
def handle_disconnect():
    logger.info('Client disconnected')

@app.route('/clear_chat', methods=['POST'])
def clear_chat():
    session_id = get_session_id()
    if session_id in chat_sessions:
        chat_sessions[session_id]['messages'] = []
    return jsonify({'message': 'Chat cleared'})

@app.route('/health')
def health_check():
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

if __name__ == '__main__':
    socketio.run(app, debug=True, host='0.0.0.0', port=5000)
