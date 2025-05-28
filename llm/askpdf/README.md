


# AskPDF - AI-Powered Document Question Answering

A powerful command-line tool that uses Amazon Bedrock's Nova models to answer questions about PDF and Markdown documents. The tool intelligently chooses between direct document submission and FAISS-based vector search depending on document size and complexity.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [AWS Setup and Credentials](#aws-setup-and-credentials)
- [Usage](#usage)
  - [Basic Usage](#basic-usage)
  - [Advanced Options](#advanced-options)
  - [Examples](#examples)
- [How It Works](#how-it-works)
  - [Direct PDF Submission](#direct-pdf-submission)
  - [FAISS Vector Search](#faiss-vector-search)
  - [Model Selection](#model-selection)
- [Configuration](#configuration)
- [Performance Optimization](#performance-optimization)
- [Troubleshooting](#troubleshooting)
  - [Installation Issues](#installation-issues)
  - [AWS Authentication Issues](#aws-authentication-issues)
  - [Runtime Errors](#runtime-errors)
  - [Performance Issues](#performance-issues)
- [Technical Details](#technical-details)
- [Contributing](#contributing)
- [License](#license)

## Features

- **Dual Processing Modes**: Automatically chooses between direct PDF submission to Nova models or FAISS-based vector search
- **Multi-Document Support**: Process single files, directories, or recursive directory searches
- **Smart Model Selection**: Automatically selects the appropriate Nova model based on document size and complexity
- **GPU Acceleration**: Supports GPU-accelerated FAISS for faster vector operations
- **Embedding Cache**: Caches embeddings to speed up repeated queries on the same documents
- **Flexible File Support**: Works with PDF and Markdown files
- **Comprehensive Error Handling**: Graceful fallbacks and detailed error messages
- **Performance Optimization**: Multiple optimization options for speed vs. accuracy tradeoffs

## Prerequisites

- Python 3.8 or higher
- AWS account with access to Amazon Bedrock
- AWS CLI configured with appropriate credentials
- NVIDIA GPU (optional, for GPU acceleration)

## Installation

### 1. Clone or Download the Script

```bash
# Download the script directly
wget https://raw.githubusercontent.com/your-repo/askpdf.py
chmod +x askpdf.py
```

### 2. Install Python Dependencies

#### Option A: Install with GPU Support (Recommended if you have NVIDIA GPU)

```bash
pip install boto3 faiss-gpu-cu12 sentence-transformers pymupdf numpy
```

#### Option B: Install with CPU-only Support

```bash
pip install boto3 faiss-cpu sentence-transformers pymupdf numpy
```

#### Option C: Install All Dependencies at Once

```bash
pip install -r requirements.txt
```

For CPU-only installation, edit `requirements.txt` and uncomment the CPU-only packages while commenting out GPU packages.

### 3. Verify Installation

```bash
python askpdf.py --help
```

## AWS Setup and Credentials

### 1. Enable Amazon Bedrock Models

1. Log into the AWS Console
2. Navigate to Amazon Bedrock
3. Go to "Model access" in the left sidebar
4. Request access to Nova models:
   - `us.amazon.nova-micro-v1:0`
   - `us.amazon.nova-lite-v1:0`
   - `us.amazon.nova-pro-v1:0`
   - `us.amazon.nova-premier-v1:0`

### 2. Configure AWS Credentials

#### Option A: AWS CLI (Recommended)

```bash
# Install AWS CLI if not already installed
pip install awscli

# Configure credentials
aws configure
```

Enter your:
- AWS Access Key ID
- AWS Secret Access Key
- Default region (e.g., `us-east-1`)
- Default output format (e.g., `json`)

#### Option B: Environment Variables

```bash
export AWS_ACCESS_KEY_ID=your_access_key_here
export AWS_SECRET_ACCESS_KEY=your_secret_key_here
export AWS_DEFAULT_REGION=us-east-1
```

#### Option C: AWS Profiles

```bash
# Configure a named profile
aws configure --profile myprofile

# Use with askpdf
python askpdf.py document.pdf "What is this about?" --profile myprofile
```

#### Option D: IAM Roles (for EC2 instances)

If running on EC2, attach an IAM role with the following permissions:
```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "bedrock:InvokeModel",
                "bedrock:InvokeModelWithResponseStream"
            ],
            "Resource": [
                "arn:aws:bedrock:*::foundation-model/us.amazon.nova-*"
            ]
        }
    ]
}
```

### 3. Test AWS Connection

```bash
# Test Bedrock access
aws bedrock list-foundation-models --region us-east-1

# Test with askpdf
python askpdf.py --help
```

## Usage

### Basic Usage

```bash
# Ask a question about a single PDF
python askpdf.py document.pdf "What are the main topics covered?"

# Ask about all PDFs in a directory
python askpdf.py /path/to/documents/ "Summarize the key findings"

# Search recursively through subdirectories
python askpdf.py /path/to/documents/ "What is the methodology?" --recursive
```

### Advanced Options

```bash
# Use specific AWS region and profile
python askpdf.py doc.pdf "Question?" --region us-west-2 --profile myprofile

# Customize model parameters
python askpdf.py doc.pdf "Question?" --temperature 0.1 --top-p 0.8 --output-tokens 2000

# Use specific model
python askpdf.py doc.pdf "Question?" --model-id us.amazon.nova-pro-v1:0

# Optimize for speed
python askpdf.py doc.pdf "Question?" --fast-embeddings --no-cache

# Force FAISS mode (disable direct PDF submission)
python askpdf.py doc.pdf "Question?" --no-faiss

# Limit number of relevant chunks
python askpdf.py doc.pdf "Question?" --top-k 5

# Custom chunk size for better granularity
python askpdf.py doc.pdf "Question?" --chunk-size 2000
```

### Examples

#### Example 1: Research Paper Analysis
```bash
python askpdf.py research_papers/ "What are the main methodologies used across these papers?" --recursive --top-k 10
```

#### Example 2: Legal Document Review
```bash
python askpdf.py contracts/ "What are the key terms and conditions?" --model-id us.amazon.nova-premier-v1:0 --temperature 0.1
```

#### Example 3: Technical Documentation
```bash
python askpdf.py technical_docs/ "How do I configure the system?" --chunk-size 4000 --fast-embeddings
```

## How It Works

### Direct PDF Submission

For smaller documents (under 25MB total), the tool can submit PDFs directly to Nova models that support document understanding:

- **Supported Models**: Nova Lite, Nova Pro, Nova Premier
- **Advantages**: Faster processing, no chunking artifacts, full document context
- **Limitations**: 25MB total size limit, PDF files only

### FAISS Vector Search

For larger documents or when direct submission isn't possible:

1. **Text Extraction**: Extracts text from PDFs using PyMuPDF
2. **Chunking**: Splits text into overlapping chunks (default 3000 characters)
3. **Embedding**: Generates vector embeddings using SentenceTransformers
4. **Indexing**: Creates FAISS index for fast similarity search
5. **Retrieval**: Finds most relevant chunks for the question
6. **Generation**: Sends relevant context to Nova model for answer generation

### Model Selection

The tool automatically selects the appropriate Nova model based on:

- Document size and estimated token count
- Model capabilities (document support)
- Context window limits
- User preferences

**Model Hierarchy**:
1. **Nova Micro** (128k tokens) - Fastest, no document support
2. **Nova Lite** (300k tokens) - Good balance, supports documents
3. **Nova Pro** (300k tokens) - Higher quality, supports documents
4. **Nova Premier** (1M tokens) - Highest quality, largest context

## Configuration

### Environment Variables

```bash
# Disable embedding cache
export ASKPDF_NO_CACHE=1

# Set custom cache directory
export ASKPDF_CACHE_DIR=/path/to/cache

# Default AWS region
export AWS_DEFAULT_REGION=us-east-1
```

### Performance Tuning

| Parameter | Default | Description | Performance Impact |
|-----------|---------|-------------|-------------------|
| `--chunk-size` | 3000 | Characters per chunk | Larger = fewer chunks, less granular |
| `--top-k` | All | Number of chunks to retrieve | Lower = faster, potentially less complete |
| `--fast-embeddings` | False | Use faster embedding model | 3x faster, slightly less accurate |
| `--no-cache` | False | Disable embedding cache | Slower repeated queries |

## Performance Optimization

### For Speed
```bash
# Use fast embeddings and limit chunks
python askpdf.py doc.pdf "Question?" --fast-embeddings --top-k 5 --chunk-size 2000
```

### For Accuracy
```bash
# Use larger chunks and premium model
python askpdf.py doc.pdf "Question?" --chunk-size 4000 --model-id us.amazon.nova-premier-v1:0
```

### For Large Documents
```bash
# Use GPU acceleration and caching
python askpdf.py large_docs/ "Question?" --recursive --chunk-size 5000
```

## Troubleshooting

### Installation Issues

#### FAISS Installation Problems

**Problem**: `ImportError: No module named 'faiss'`

**Solutions**:
```bash
# For GPU systems
pip uninstall faiss-gpu-cu12 faiss-cpu
pip install faiss-gpu-cu12

# For CPU-only systems
pip uninstall faiss-gpu-cu12 faiss-cpu
pip install faiss-cpu

# If still failing, try conda
conda install -c conda-forge faiss-gpu
# or
conda install -c conda-forge faiss-cpu
```

#### PyMuPDF Installation Issues

**Problem**: `ImportError: No module named 'fitz'`

**Solution**:
```bash
pip uninstall pymupdf
pip install pymupdf==1.23.26
```

#### SentenceTransformers Issues

**Problem**: Model download failures or CUDA errors

**Solutions**:
```bash
# Clear cache and reinstall
rm -rf ~/.cache/torch/sentence_transformers/
pip uninstall sentence-transformers
pip install sentence-transformers

# For CUDA issues
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

### AWS Authentication Issues

#### Problem: `NoCredentialsError`

**Solutions**:
1. Check AWS credentials:
   ```bash
   aws sts get-caller-identity
   ```

2. Reconfigure AWS CLI:
   ```bash
   aws configure
   ```

3. Check environment variables:
   ```bash
   echo $AWS_ACCESS_KEY_ID
   echo $AWS_SECRET_ACCESS_KEY
   ```

#### Problem: `AccessDeniedException`

**Solutions**:
1. Verify Bedrock model access in AWS Console
2. Check IAM permissions
3. Ensure you're using the correct region:
   ```bash
   python askpdf.py doc.pdf "Question?" --region us-east-1
   ```

#### Problem: Model not available in region

**Solution**:
```bash
# Try different regions where Nova models are available
python askpdf.py doc.pdf "Question?" --region us-east-1
python askpdf.py doc.pdf "Question?" --region us-west-2
```

### Runtime Errors

#### Problem: `ValidationException: Input is too long`

**Solutions**:
1. Use smaller chunk size:
   ```bash
   python askpdf.py doc.pdf "Question?" --chunk-size 2000 --top-k 5
   ```

2. Force FAISS mode:
   ```bash
   python askpdf.py doc.pdf "Question?" --no-faiss
   ```

3. Use larger model:
   ```bash
   python askpdf.py doc.pdf "Question?" --model-id us.amazon.nova-premier-v1:0
   ```

#### Problem: `ReadTimeoutError`

**Solutions**:
1. Reduce document size or complexity
2. Use direct PDF submission for smaller files:
   ```bash
   python askpdf.py small_doc.pdf "Question?" --no-faiss
   ```

3. Break large documents into smaller files

#### Problem: Memory errors with FAISS

**Solutions**:
1. Use CPU-only FAISS:
   ```bash
   pip uninstall faiss-gpu-cu12
   pip install faiss-cpu
   ```

2. Reduce chunk size:
   ```bash
   python askpdf.py doc.pdf "Question?" --chunk-size 1500
   ```

3. Process fewer documents at once

### Performance Issues

#### Problem: Slow embedding generation

**Solutions**:
1. Use fast embeddings:
   ```bash
   python askpdf.py doc.pdf "Question?" --fast-embeddings
   ```

2. Enable GPU acceleration (if available)

3. Use embedding cache (enabled by default)

#### Problem: Slow vector search

**Solutions**:
1. Limit number of chunks:
   ```bash
   python askpdf.py doc.pdf "Question?" --top-k 10
   ```

2. Use GPU-accelerated FAISS

3. Increase chunk size to reduce total chunks:
   ```bash
   python askpdf.py doc.pdf "Question?" --chunk-size 5000
   ```

#### Problem: High memory usage

**Solutions**:
1. Process documents individually instead of in batches
2. Use smaller chunk sizes
3. Clear embedding cache:
   ```bash
   rm -rf .embedding_cache/
   ```

### Common Error Messages and Solutions

| Error Message | Cause | Solution |
|---------------|-------|----------|
| `No PDF files found` | Wrong path or no PDFs | Check file path and extensions |
| `FAISS installation appears corrupted` | Broken FAISS install | Reinstall FAISS |
| `Model access denied` | No Bedrock model access | Enable models in AWS Console |
| `Region not supported` | Model not in region | Use `--region us-east-1` |
| `Input is too long` | Document too large | Use `--chunk-size 2000 --top-k 5` |
| `Connection timeout` | Network/AWS issues | Check internet and AWS status |

## Technical Details

### Supported File Types
- **PDF**: All standard PDF files (text-based)
- **Markdown**: `.md` files with UTF-8 encoding

### Model Limits
| Model | Context Window | Document Support | Best For |
|-------|----------------|------------------|----------|
| Nova Micro | 128k tokens | No | Small text queries |
| Nova Lite | 300k tokens | Yes | General document QA |
| Nova Pro | 300k tokens | Yes | Complex analysis |
| Nova Premier | 1M tokens | Yes | Large documents |

### Performance Benchmarks
- **Small PDFs** (< 5MB): 10-30 seconds
- **Medium PDFs** (5-20MB): 30-120 seconds  
- **Large document sets**: 2-10 minutes
- **GPU acceleration**: 2-3x faster embedding generation

### Cache Behavior
- Embeddings cached by file hash and chunk parameters
- Cache location: `.embedding_cache/` directory
- Cache invalidated when file content or chunk size changes
- Disable with `--no-cache` flag

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is licensed under the MIT License. See LICENSE file for details.

---

**Need Help?** 
- Check the troubleshooting section above
- Review AWS Bedrock documentation
- Open an issue on GitHub with detailed error messages and system information
