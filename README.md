# dptests

Just a bit of test code here and there


# AI Chatbot with Document Upload

A modern Flask-based chatbot application that integrates with AWS Bedrock Nova Lite model and supports document uploads for context-aware conversations.

## Features

- 🤖 **AI-Powered Chat**: Uses AWS Bedrock Nova Lite v1 model
- 📁 **Document Upload**: Support for multiple file types (TXT, PDF, DOC, DOCX, MD, JSON, CSV)
- 💬 **Real-time Communication**: WebSocket-based chat using Socket.IO
- 📱 **Responsive Design**: Modern Web 2.0 interface that works on all devices
- 🔒 **Secure File Handling**: File validation and secure storage
- 💾 **Session Management**: Maintains chat history and document context

## Prerequisites

- Python 3.8+
- AWS Account with Bedrock access
- AWS credentials configured

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd chatbot-app
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Set up environment variables:
```bash
cp .env.example .env
# Edit .env with your AWS credentials
```

4. Create uploads directory:
```bash
mkdir uploads
```

## Configuration

### AWS Setup

1. Ensure you have AWS Bedrock access in your region
2. Configure your AWS credentials either through:
   - Environment variables (AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY)
   - AWS CLI (`aws configure`)
   - IAM roles (if running on EC2)

### Required AWS Permissions

Your AWS user/role needs the following permissions:
```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "bedrock:InvokeModel"
            ],
            "Resource": "arn:aws:bedrock:*:*:model/us.amazon.nova-lite-v1:0"
        }
    ]
}
```

## Usage

1. Start the application:
```bash
python app.py
```

2. Open your browser and navigate to `http://localhost:5000`

3. Upload documents (optional) by dragging and dropping or clicking the upload area

4. Start chatting with the AI assistant

## File Upload Limits

- Maximum file size: 1GB total
- Supported formats: TXT, PDF, DOC, DOCX, MD, JSON, CSV
- Multiple files can be uploaded simultaneously

## API Endpoints

- `GET /` - Main chat interface
- `POST /upload` - File upload endpoint
- `POST /clear_chat` - Clear chat history
- `GET /health` - Health check endpoint

## WebSocket Events

- `send_message` - Send a message to the AI
- `message_response` - Receive messages from the AI
- `connect` - Client connection established
- `disconnect` - Client disconnected

## Production Deployment

For production deployment, consider:

1. Use a production WSGI server like Gunicorn:
```bash
pip install gunicorn
gunicorn --worker-class eventlet -w 1 app:app
```

2. Set up a reverse proxy (nginx)
3. Use a proper database for session storage
4. Configure proper logging
5. Set up SSL/TLS certificates
6. Use environment-specific configuration

## Security Considerations

- File uploads are validated and stored securely
- Session management uses secure cookies
- AWS credentials should be properly secured
- Consider implementing rate limiting for production use

## Troubleshooting

### Common Issues

1. **AWS Bedrock Access Denied**: Ensure your AWS credentials have the required permissions
2. **File Upload Errors**: Check file size limits and supported formats
3. **Connection Issues**: Verify WebSocket connectivity and firewall settings

### Logs

The application logs important events. Check the console output for debugging information.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License.
# BedBot - AI Chat Assistant

BedBot is a powerful AI chat assistant powered by AWS Bedrock Nova that enables context-aware conversations with document support. Upload documents (PDF, TXT, DOC, DOCX, MD, JSON, CSV) and chat with an AI that understands their content.

## Table of Contents

- [Features](#features)
- [Use Cases](#use-cases)
- [Prerequisites](#prerequisites)
- [AWS Setup](#aws-setup)
  - [Creating an AWS Profile](#creating-an-aws-profile)
  - [Requesting Bedrock Model Access](#requesting-bedrock-model-access)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Basic Chat](#basic-chat)
  - [Document Upload](#document-upload)
  - [Voice Input](#voice-input)
  - [Command Line Options](#command-line-options)
- [Architecture](#architecture)
  - [Storage Modes](#storage-modes)
  - [Document Processing](#document-processing)
- [Troubleshooting](#troubleshooting)
- [Security Considerations](#security-considerations)
- [Contributing](#contributing)

## Features

- **Context-Aware Chat**: Upload documents and chat with AI that understands their content
- **Multi-Format Support**: PDF, TXT, DOC, DOCX, MD, JSON, CSV files
- **Dual Storage Modes**: Local filesystem or AWS S3 bucket storage
- **Voice Input**: Speech-to-text functionality for hands-free interaction
- **Document Comparison**: Upload multiple documents for comparative analysis
- **PDF Merging**: Automatically merges multiple PDFs in S3 mode for efficient processing
- **Real-time Chat**: Modern web interface with typing indicators and markdown support
- **Session Management**: Persistent sessions with file management
- **Automatic Cleanup**: Temporary files and S3 buckets are cleaned up automatically

## Use Cases

### Document Analysis
- **Legal Document Review**: Upload contracts, agreements, or legal documents for analysis and comparison
- **Research Paper Analysis**: Compare multiple research papers, extract key findings, and identify relationships
- **Technical Documentation**: Get explanations of complex technical documents, manuals, or specifications
- **Financial Report Analysis**: Upload financial statements, reports, or proposals for detailed analysis

### Content Creation
- **Report Writing**: Upload source materials and get help writing comprehensive reports
- **Content Summarization**: Extract key points and summaries from lengthy documents
- **Comparative Analysis**: Upload multiple documents to identify differences, similarities, and recommendations

### Educational Support
- **Study Aid**: Upload textbooks, papers, or course materials for Q&A sessions
- **Research Assistance**: Get help understanding complex academic papers or research materials
- **Document Synthesis**: Combine information from multiple sources into coherent explanations

### Business Intelligence
- **Proposal Evaluation**: Upload RFPs and proposals for detailed comparison and evaluation
- **Policy Analysis**: Compare policy documents, regulations, or compliance materials
- **Market Research**: Analyze market reports, competitor analysis, and industry documents

## Prerequisites

- Python 3.8 or higher
- AWS Account with Bedrock access
- Modern web browser (Chrome, Firefox, Safari, Edge)

## AWS Setup

### Creating an AWS Profile

If you haven't set up AWS yet, you need AWS credentials from your AWS administrator. Ideally, ask for a dedicated IAM user with `AmazonBedrockFullAccess` permissions.

#### Install AWS CLI

If the `aws` command is not available, install it:

**Linux/Mac:**
```bash
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
```

**Windows:**
Download and run the AWS CLI MSI installer from the [official AWS documentation](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

#### Configure AWS Profile

Create a new AWS profile specifically for Bedrock:

```bash
aws configure --profile bedbot

AWS Access Key ID [None]: YOUR_ACCESS_KEY
AWS Secret Access Key [None]: YOUR_SECRET_KEY
Default region name [None]: us-west-2
Default output format [None]: json
```

#### Enterprise SSO Setup

If your organization uses AWS SSO, configure it instead:

```bash
aws configure sso --profile bedbot

SSO session name [None]: my-sso
SSO start URL [None]: https://your-org.awsapps.com/start
SSO region [None]: us-west-2
SSO registration scopes [None]: sso:account:access
```

### Requesting Bedrock Model Access

**Important**: In addition to Bedrock permissions, you need to request access to specific LLM models in your AWS account. This is a one-time setup per AWS account.

1. Go to the AWS Bedrock console
2. Navigate to "Model access" in the left sidebar
3. Click "Request model access"
4. Select the models you want to use (recommended: Nova Premier, Nova Lite, Claude models)
5. Submit the request

[Watch this 2-minute video on requesting model access](https://www.youtube.com/watch?v=WWHo7Awy0sQ)

## Installation

1. **Clone or download BedBot:**
```bash
git clone <repository-url>
cd bedbot
```

2. **Install Python dependencies:**
```bash
pip install -r requirements.txt
```

3. **Set environment variables (optional):**
```bash
export AWS_PROFILE=bedbot
export SECRET_KEY=your-secret-key-here
```

## Configuration

### Environment Variables

- `AWS_PROFILE`: AWS profile to use (default: default profile)
- `SECRET_KEY`: Flask session secret key (auto-generated if not set)
- `WERKZEUG_RUN_MAIN`: Internal Flask variable (don't modify)

### Model Selection

BedBot supports various AWS Bedrock models:
- `us.amazon.nova-premier-v1:0` (default) - Best balance of capability and speed
- `us.amazon.nova-lite-v1:0` - Faster, more economical
- `us.anthropic.claude-3-5-sonnet-20241022-v2:0` - High capability
- `us.anthropic.claude-3-5-haiku-20241022-v1:0` - Fast and efficient

## Usage

### Basic Chat

Start BedBot with default settings:

```bash
python bedbot.py
```

Open your browser to `http://localhost:5000` and start chatting!

### Document Upload

1. **Drag and drop** files onto the upload area, or **click** to select files
2. Supported formats: TXT, PDF, DOC, DOCX, MD, JSON, CSV
3. Maximum file size: 4.5MB per file (local mode), unlimited (S3 mode)
4. Maximum files: 1000 per session

### Voice Input

1. Click the **microphone button** in the chat interface
2. Allow microphone access when prompted
3. Speak your message
4. Click the **stop button** or the microphone button again to finish
5. Your speech will be converted to text automatically

### Command Line Options

```bash
# Use local filesystem instead of S3
python bedbot.py --no-bucket

# Enable debug mode for detailed logging
python bedbot.py --debug

# Use a specific Bedrock model
python bedbot.py --model us.anthropic.claude-3-5-sonnet-20241022-v2:0

# Combine options
python bedbot.py --no-bucket --debug --model us.amazon.nova-lite-v1:0
```

## Architecture

### Storage Modes

#### S3 Mode (Default)
- Creates temporary S3 bucket for file storage
- Supports unlimited file sizes
- Automatically merges multiple PDFs for efficient processing
- Bucket is automatically deleted on shutdown
- Better for large documents and production use

#### Local Mode (`--no-bucket`)
- Stores files in temporary local directories
- Files larger than 4.5MB are truncated
- Suitable for development and smaller documents
- All temporary files cleaned up on shutdown

### Document Processing

1. **Text Files**: Content is extracted and added to conversation context
2. **PDF Files**: Processed using AWS Bedrock Nova's native document understanding
3. **Multiple PDFs**: Automatically merged in S3 mode for comparative analysis
4. **Large Files**: Truncated in local mode, preserved in S3 mode

## Troubleshooting

### Common Issues

#### "Bedrock client not initialized"
**Problem**: AWS credentials not configured properly.

**Solutions**:
- Verify AWS profile: `aws sts get-caller-identity --profile bedbot`
- Check credentials: `aws configure list --profile bedbot`
- Ensure Bedrock permissions are granted to your IAM user/role

#### "Failed to create S3 bucket"
**Problem**: S3 bucket creation failed, falling back to local mode.

**Solutions**:
- Check S3 permissions in your IAM policy
- Verify AWS region supports S3 and Bedrock
- Use `--no-bucket` flag to force local mode
- Check AWS service quotas for S3 buckets

#### "Model access denied"
**Problem**: Requested Bedrock model not accessible.

**Solutions**:
- Request model access in AWS Bedrock console
- Wait for model access approval (can take a few minutes)
- Try a different model: `--model us.amazon.nova-lite-v1:0`
- Check model availability in your AWS region

#### "Upload failed: File too large"
**Problem**: File exceeds size limits.

**Solutions**:
- Use S3 mode (remove `--no-bucket` flag) for larger files
- Split large documents into smaller sections
- Use PDF compression tools to reduce file size

#### Voice input not working
**Problem**: Speech recognition not functioning.

**Solutions**:
- Use Chrome or Edge browser (better speech recognition support)
- Allow microphone access when prompted
- Check browser microphone permissions in settings
- Ensure microphone is working in other applications

### Debug Mode

Enable debug mode for detailed logging:

```bash
python bedbot.py --debug
```

Debug mode shows:
- AWS API request/response details
- File processing information
- S3 bucket operations
- Session management details

### Log Analysis

Check the console output for:
- AWS authentication status
- File upload progress
- Model response times
- Error messages with stack traces

### Performance Issues

#### Slow responses
- Try a faster model: `--model us.amazon.nova-lite-v1:0`
- Reduce document size or number of documents
- Check AWS region latency
- Verify internet connection stability

#### High AWS costs
- Use Nova Lite instead of Nova Premier
- Limit document sizes and quantities
- Monitor AWS Bedrock usage in AWS console
- Consider using local mode for development

## Security Considerations

### Data Privacy
- **Temporary Storage**: All files are stored temporarily and automatically deleted
- **Session Isolation**: Each session has isolated file storage
- **No Persistent Data**: No chat history or files are permanently stored
- **S3 Encryption**: S3 buckets use server-side encryption (AES256)

### Access Control
- **IAM Permissions**: Use least-privilege IAM policies
- **Network Security**: Run on localhost by default
- **Session Security**: Flask sessions are server-side and encrypted
- **File Validation**: Only allowed file types are processed

### Best Practices
- Use dedicated AWS IAM user for BedBot
- Regularly rotate AWS access keys
- Monitor AWS CloudTrail for API usage
- Don't upload sensitive documents in shared environments
- Use HTTPS in production deployments

### Production Deployment
For production use:
- Set strong `SECRET_KEY` environment variable
- Use HTTPS with proper SSL certificates
- Implement authentication/authorization
- Configure proper logging and monitoring
- Use AWS VPC for network isolation

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make your changes and test thoroughly
4. Commit with descriptive messages: `git commit -m "Add feature description"`
5. Push to your fork: `git push origin feature-name`
6. Create a Pull Request

### Development Setup

```bash
# Clone repository
git clone <repository-url>
cd bedbot

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate     # Windows

# Install dependencies
pip install -r requirements.txt

# Run in debug mode
python bedbot.py --debug --no-bucket
```

### Testing

Test with various file types and sizes:
```bash
# Test with sample documents
python bedbot.py --debug --no-bucket

# Upload test files through the web interface
# Try different models and configurations
```

---

**BedBot** - Empowering conversations with AI and documents
