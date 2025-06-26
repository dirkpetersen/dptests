

# BedBot - AI Chat Assistant

BedBot is a powerful AI chat assistant powered by AWS Bedrock Nova that allows you to have intelligent conversations and analyze documents. It supports both text-based conversations and document analysis with PDF processing capabilities.


![image](https://github.com/user-attachments/assets/b9894fb3-2204-4a9f-a165-5de122f16bb1)


## Table of Contents

- [BedBot - AI Chat Assistant](#bedbot---ai-chat-assistant)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Use Cases](#use-cases)
    - [Document Analysis and Comparison](#document-analysis-and-comparison)
    - [Content Creation and Review](#content-creation-and-review)
    - [Educational and Training](#educational-and-training)
    - [Business Intelligence](#business-intelligence)
  - [Prerequisites](#prerequisites)
  - [AWS Setup](#aws-setup)
    - [Creating a New AWS Profile](#creating-a-new-aws-profile)
    - [Requesting Bedrock Model Access](#requesting-bedrock-model-access)
  - [Installation](#installation)
  - [Configuration](#configuration)
    - [Environment Variables](#environment-variables)
    - [Model Selection](#model-selection)
  - [Usage](#usage)
    - [Basic Startup](#basic-startup)
    - [Basic Chat](#basic-chat)
    - [Document Upload and Analysis](#document-upload-and-analysis)
    - [Voice Input](#voice-input)
  - [Command Line Options](#command-line-options)
  - [Storage Modes](#storage-modes)
    - [Local Filesystem Mode (`--no-bucket`)](#local-filesystem-mode---no-bucket)
    - [S3 Storage Mode (default)](#s3-storage-mode-default)
  - [Troubleshooting](#troubleshooting)
    - [Common Issues](#common-issues)
      - [1. "Bedrock client not initialized"](#1-bedrock-client-not-initialized)
      - [2. "Model access denied" or "ValidationException"](#2-model-access-denied-or-validationexception)
      - [3. "Failed to create S3 bucket"](#3-failed-to-create-s3-bucket)
      - [4. File upload fails](#4-file-upload-fails)
      - [5. Voice input not working](#5-voice-input-not-working)
    - [Debug Mode](#debug-mode)
    - [Log Analysis](#log-analysis)
  - [Advanced Features](#advanced-features)
    - [Session Management](#session-management)
    - [Document Context](#document-context)
    - [PDF Processing](#pdf-processing)
    - [Security Features](#security-features)
  - [Security Considerations](#security-considerations)
    - [Data Privacy](#data-privacy)
    - [AWS Permissions](#aws-permissions)
    - [Network Security](#network-security)
    - [Recommended IAM Policy](#recommended-iam-policy)
  - [Contributing](#contributing)
    - [Development Setup](#development-setup)

## Features

- **AI-Powered Conversations**: Chat with AWS Bedrock Nova models for intelligent responses
- **Document Analysis**: Upload and analyze multiple document types (PDF, TXT, DOC, DOCX, MD, JSON, CSV)
- **PDF Processing**: Advanced PDF handling with automatic merging and content extraction
- **Voice Input**: Speech-to-text functionality for hands-free interaction
- **Dual Storage Modes**: Choose between local filesystem or AWS S3 for document storage
- **Session Management**: Persistent chat history and document context within sessions
- **Modern Web Interface**: Responsive, mobile-friendly chat interface
- **Real-time Processing**: Live typing indicators and instant responses
- **File Management**: Easy upload, preview, and removal of documents

## Use Cases

### Document Analysis and Comparison
- **Legal Document Review**: Upload contracts, agreements, or legal documents for analysis and comparison
- **Research Paper Analysis**: Compare multiple research papers, extract key findings, and identify relationships
- **Technical Documentation**: Analyze technical specifications, requirements documents, and implementation guides
- **Financial Report Analysis**: Review financial statements, compare quarterly reports, and extract insights

### Content Creation and Review
- **Content Summarization**: Upload long documents for concise summaries and key point extraction
- **Writing Assistance**: Get help with writing, editing, and improving documents
- **Translation and Localization**: Analyze documents in different languages and get translation assistance
- **Compliance Checking**: Compare documents against regulatory requirements or standards

### Educational and Training
- **Study Material Analysis**: Upload textbooks, papers, and study materials for Q&A sessions
- **Training Document Review**: Analyze training materials and create interactive learning experiences
- **Curriculum Development**: Compare educational standards and develop comprehensive curricula

### Business Intelligence
- **Market Research**: Analyze market reports, competitor analysis, and industry trends
- **Policy Analysis**: Review company policies, procedures, and compliance documents
- **Strategic Planning**: Analyze business plans, strategic documents, and performance reports

## Prerequisites

- Python 3.8 or higher
- AWS Account with Bedrock access
- AWS CLI installed and configured
- Modern web browser with JavaScript enabled

## AWS Setup

### Creating a New AWS Profile

If you haven't set up AWS yet, you'll need AWS credentials from your AWS administrator. Ideally, ask for a dedicated IAM user with `AmazonBedrockFullAccess` permissions.

1. **Install AWS CLI** (if not already installed):
   
   **Linux/Mac:**
   ```bash
   curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
   unzip awscliv2.zip
   sudo ./aws/install
   ```
   
   **Windows:**
   ```powershell
   msiexec.exe /i https://awscli.amazonaws.com/AWSCLIV2.msi
   ```

2. **Configure AWS Profile:**
   ```bash
   aws configure --profile bedbot
   ```
   
   Enter your credentials:
   ```
   AWS Access Key ID [None]: YOUR_ACCESS_KEY
   AWS Secret Access Key [None]: YOUR_SECRET_KEY
   Default region name [None]: us-east-1
   Default output format [None]: json
   ```

3. **For Enterprise/SSO Users:**
   If your organization uses AWS SSO, configure it instead:
   ```bash
   aws configure sso --profile bedbot
   ```

### Requesting Bedrock Model Access

**Important**: In addition to Bedrock permissions, you need to request access to specific models in your AWS account.

1. Go to the AWS Bedrock console
2. Navigate to "Model access" in the left sidebar
3. Click "Request model access"
4. Select the models you want to use (recommended: Nova Premier, Nova Lite, Claude models)
5. Submit the request (usually approved within minutes)

**Video Guide**: [2-minute video on requesting model access](https://www.youtube.com/watch?v=WWHo7Awy0sQ)

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

3. **Verify AWS configuration:**
   ```bash
   aws bedrock list-foundation-models --profile bedbot --region us-east-1
   ```

## Configuration

### Environment Variables

Create a `.env` file in the project directory (optional):

```bash
AWS_PROFILE=bedbot
AWS_DEFAULT_REGION=us-east-1
SECRET_KEY=your-secret-key-change-this
```

### Model Selection

BedBot supports various Bedrock models. You can specify the model using the `--model` parameter:

- `us.amazon.nova-premier-v1:0` (default) - Most capable Nova model
- `us.amazon.nova-lite-v1:0` - Faster, more cost-effective
- `us.anthropic.claude-3-5-sonnet-20241022-v2:0` - Anthropic's Claude
- `us.anthropic.claude-3-5-haiku-20241022-v1:0` - Fast Claude model

## Usage

### Basic Startup

**Local filesystem mode (default):**
```bash
python bedbot.py
```

**S3 storage mode:**
```bash
python bedbot.py --no-bucket=false
```

**With custom model:**
```bash
python bedbot.py --model us.amazon.nova-lite-v1:0
```

**Debug mode:**
```bash
python bedbot.py --debug
```

### Basic Chat

1. Open your browser to `http://localhost:5000`
2. Type your message in the chat input
3. Press Enter or click the send button
4. BedBot will respond using the configured Bedrock model

### Document Upload and Analysis

1. **Upload Documents:**
   - Click the upload area or drag and drop files
   - Supported formats: PDF, TXT, DOC, DOCX, MD, JSON, CSV
   - Maximum file size: 4.5MB per file (local mode) or unlimited (S3 mode)
   - Maximum files: 1000 per session

2. **Document Processing:**
   - PDFs are processed using Nova's native document understanding
   - Text files are extracted and added to conversation context
   - Multiple PDFs are automatically merged in S3 mode for better analysis

3. **Ask Questions:**
   ```
   "What are the key points in the uploaded document?"
   "Compare the requirements in document A with the proposal in document B"
   "Summarize the main findings from all uploaded documents"
   ```

### Voice Input

1. Click the microphone button
2. Speak your message clearly
3. Click the stop button or the microphone again to finish
4. Your speech will be converted to text and can be edited before sending

**Note**: Voice input requires a modern browser with Web Speech API support (Chrome, Edge, Safari).

## Command Line Options

```bash
python bedbot.py [OPTIONS]
```

**Options:**
- `--no-bucket`: Use local filesystem instead of S3 bucket (default: false)
- `--debug`: Enable debug mode to print API messages and detailed logs
- `--model MODEL`: Specify Bedrock model to use (default: us.amazon.nova-premier-v1:0)

**Examples:**
```bash
# Local mode with debug
python bedbot.py --no-bucket --debug

# S3 mode with Claude model
python bedbot.py --model us.anthropic.claude-3-5-sonnet-20241022-v2:0

# Fast Nova Lite model
python bedbot.py --model us.amazon.nova-lite-v1:0 --debug
```

## Storage Modes

### Local Filesystem Mode (`--no-bucket`)
- **Pros**: No AWS S3 costs, faster for small files, works offline for chat
- **Cons**: Limited file size (4.5MB), no persistence across server restarts
- **Use Case**: Development, testing, small document analysis

### S3 Storage Mode (default)
- **Pros**: Unlimited file sizes, automatic cleanup, better for large documents
- **Cons**: Requires S3 permissions, incurs AWS storage costs
- **Use Case**: Production use, large document analysis, team collaboration

## Troubleshooting

### Common Issues

#### 1. "Bedrock client not initialized"
**Cause**: AWS credentials not configured properly
**Solution**:
```bash
# Check AWS configuration
aws configure list --profile bedbot

# Test Bedrock access
aws bedrock list-foundation-models --profile bedbot --region us-east-1
```

#### 2. "Model access denied" or "ValidationException"
**Cause**: Model access not requested in Bedrock console
**Solution**:
1. Go to AWS Bedrock console → Model access
2. Request access to the model you're trying to use
3. Wait for approval (usually immediate)

#### 3. "Failed to create S3 bucket"
**Cause**: Insufficient S3 permissions or bucket name conflicts
**Solution**:
```bash
# Use local mode instead
python bedbot.py --no-bucket

# Or check S3 permissions
aws s3 ls --profile bedbot
```

#### 4. File upload fails
**Cause**: File too large or unsupported format
**Solution**:
- Check file size (4.5MB limit in local mode)
- Verify file format is supported
- Try S3 mode for larger files: remove `--no-bucket` flag

#### 5. Voice input not working
**Cause**: Browser doesn't support Web Speech API or microphone permissions denied
**Solution**:
- Use Chrome, Edge, or Safari
- Allow microphone access when prompted
- Check browser console for errors

### Debug Mode

Enable debug mode for detailed logging:
```bash
python bedbot.py --debug
```

This will show:
- AWS API calls and responses
- File processing details
- Session management information
- Error stack traces

### Log Analysis

Check the console output for:
- AWS credential issues
- Model access problems
- File processing errors
- Network connectivity issues

## Advanced Features

### Session Management
- Each browser session maintains separate chat history and document context
- Sessions persist until explicitly cleared or server restart (local mode)
- S3 mode provides better session persistence

### Document Context
- Uploaded documents become part of the conversation context
- BedBot can reference and compare multiple documents
- Context is maintained throughout the session

### PDF Processing
- Automatic PDF merging for multiple uploads (S3 mode)
- Native PDF understanding using Nova's document capabilities
- Intelligent content extraction and analysis

### Security Features
- Server-side session management
- Automatic cleanup of temporary files
- S3 bucket encryption and access controls
- No persistent storage of sensitive data

## Security Considerations

### Data Privacy
- Documents are processed temporarily and cleaned up automatically
- S3 buckets are created with private access and encryption
- No data is permanently stored unless explicitly configured

### AWS Permissions
- Use dedicated IAM users with minimal required permissions
- Regularly rotate AWS access keys
- Monitor AWS CloudTrail for API usage

### Network Security
- BedBot runs on localhost by default
- Use HTTPS in production deployments
- Consider VPN or private network access for sensitive documents

### Recommended IAM Policy

For the IAM user account running BedBot create 2 new inline policies 

Policy: UseBedrock

```json
{
	"Version": "2012-10-17",
	"Statement": [
		{
			"Effect": "Allow",
			"Action": [
				"bedrock:InvokeModel",
				"bedrock:InvokeModelWithResponseStream",
				"bedrock:ListFoundationModels"
			],
			"Resource": "*"
		}
	]
}
```

Policy: AllowS3BucketsWithPrefix_bedbot

```json
{
	"Version": "2012-10-17",
	"Statement": [
		{
			"Effect": "Allow",
			"Action": "s3:*",
			"Resource": [
				"arn:aws:s3:::bedbot-*",
				"arn:aws:s3:::bedbot-*/*"
			]
		}
	]
}
```


## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

### Development Setup
```bash
# Install development dependencies
pip install -r requirements.txt

# Run in debug mode
python bedbot.py --debug --no-bucket

# Run tests (if available)
python -m pytest tests/
```

---

**Need Help?** 
- Check the [Troubleshooting](#troubleshooting) section
- Enable `--debug` mode for detailed logs
- Review AWS Bedrock documentation
- Check AWS service status if experiencing connectivity issues

**Version**: 1.0.0  
**License**: MIT  
**AWS Services**: Bedrock, S3 (optional)
