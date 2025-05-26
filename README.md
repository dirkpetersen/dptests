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
