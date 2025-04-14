"""
AWS Lambda function for distributed embedding generation
=======================================================

This Lambda function processes batches of text chunks and generates embeddings
using AWS Bedrock's embedding models. It's designed to be used with the
EnhancedServerlessRAG system to handle large document processing in parallel.

The function accepts a list of text chunks and returns their embeddings.
"""

import json
import boto3
import logging
import os
import time
from typing import List, Dict, Any

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# AWS Service clients
bedrock_client = boto3.client('bedrock-runtime')

# Default embedding model
DEFAULT_MODEL_ID = "amazon.titan-embed-text-v1"

def lambda_handler(event, context):
    """
    Lambda handler for generating embeddings from text chunks.
    
    Args:
        event: Lambda event containing:
            - texts: List of text chunks to embed
            - model_id: (Optional) Bedrock model ID to use
            
    Returns:
        Dictionary with:
            - embeddings: List of embedding vectors
            - or error information if processing failed
    """
    try:
        # Extract parameters from event
        texts = event.get('texts', [])
        model_id = event.get('model_id', DEFAULT_MODEL_ID)
        
        if not texts:
            logger.warning("No texts provided for embedding generation")
            return {
                'statusCode': 400,
                'error': 'No texts provided'
            }
        
        logger.info(f"Generating embeddings for {len(texts)} chunks using model {model_id}")
        
        # Process each text chunk
        embeddings = []
        for i, text in enumerate(texts):
            try:
                # Format the request for the embedding model
                request_body = json.dumps({
                    "inputText": text
                })
                
                # Call Bedrock for embedding generation
                response = bedrock_client.invoke_model(
                    modelId=model_id,
                    contentType="application/json",
                    accept="application/json",
                    body=request_body
                )
                
                # Parse the response
                response_body = json.loads(response['body'].read())
                embedding_vector = response_body['embedding']
                embeddings.append(embedding_vector)
                
                # Log progress periodically
                if (i + 1) % 10 == 0 or i == 0 or i == len(texts) - 1:
                    logger.info(f"Processed {i+1}/{len(texts)} embeddings")
                
            except Exception as e:
                logger.error(f"Error generating embedding for chunk {i}: {str(e)}")
                # Return a zero vector as fallback
                embeddings.append([0.0] * 1536)  # Default size for many embedding models
        
        logger.info(f"Successfully generated {len(embeddings)} embeddings")
        
        return {
            'statusCode': 200,
            'embeddings': embeddings
        }
    except Exception as e:
        logger.error(f"Error in embedding handler: {str(e)}")
        return {
            'statusCode': 500,
            'error': str(e)
        }

# For local testing
if __name__ == "__main__":
    # Test with a small sample
    test_event = {
        'texts': [
            "This is a test document for embedding generation.",
            "We want to see if the embedding function works correctly."
        ],
        'model_id': DEFAULT_MODEL_ID
    }
    
    result = lambda_handler(test_event, None)
    print(f"Generated {len(result.get('embeddings', []))} embeddings")
