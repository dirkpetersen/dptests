import os

# AWS Bedrock Configuration
BEDROCK_REGION = "us-east-1"
MODEL_ID = "anthropic.claude-3-sonnet-20240229-v1:0"  # Updated model ID
ANTHROPIC_VERSION = "bedrock-2023-05-31"
ASSUMED_ROLE = os.getenv('AWS_ROLE_ARN')  # Role ARN for Bedrock access
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
