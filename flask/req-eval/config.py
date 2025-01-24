import os

# AWS Bedrock Configuration
BEDROCK_REGION = "us-west-2"  # Default region if not configured
MODEL_ID = "anthropic.claude-3-sonnet-20240229-v1:0"
ANTHROPIC_VERSION = "bedrock-2023-05-31"
MAX_RETRIES = 10

# File Configuration
MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB

# Text Processing Configuration
CHUNK_SIZE = 1000
CHUNK_OVERLAP = 100
SIMILAR_CHUNKS_COUNT = 2

# Response Configuration
MAX_TOKENS = 1000
TEMPERATURE = 0
