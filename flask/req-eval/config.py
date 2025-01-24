import os
import logging

# Logging Configuration
LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# AWS Bedrock Configuration
BEDROCK_REGION = "us-west-2"  # Default region if not configured
MODEL_ID = "anthropic.claude-3-5-haiku-20241022-v1:0"  # Using Haiku model which accepts shorter input
ANTHROPIC_VERSION = "bedrock-2023-05-31"
MAX_RETRIES = 10

# File Configuration
MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB

# Response Configuration
MAX_TOKENS = 1000
TEMPERATURE = 0
MAX_CHARS_PER_DOC = 6000  # Reduced from 12000 to work with Haiku model

# Cache Configuration
CACHE_DIR = "policy_cache"  # Directory for policy caches
COOKIE_NAME = "user_id"  # Cookie name for user identification
COOKIE_MAX_AGE = 30 * 24 * 60 * 60  # 30 days in seconds
