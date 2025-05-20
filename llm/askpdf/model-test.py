#!/usr/bin/env python3

import argparse
import boto3
import json
import sys

# Default AWS Region
AWS_REGION = "us-east-1"

# Default model ID
DEFAULT_MODEL_ID = "us.amazon.nova-pro-v1:0"

# Default inference parameters
DEFAULT_MAX_TOKENS = 1024
DEFAULT_TOP_P = 0.9
DEFAULT_TEMPERATURE = 0.2

def main():
    parser = argparse.ArgumentParser(
        description="Test an Amazon Bedrock model with a simple prompt."
    )
    parser.add_argument(
        "--model-id",
        type=str,
        default=DEFAULT_MODEL_ID,
        help=f"Bedrock model ID to test. Default: {DEFAULT_MODEL_ID}"
    )
    parser.add_argument(
        "--prompt",
        type=str,
        default="Explain what makes a good prompt for an LLM in 3 bullet points.",
        help="The prompt to send to the model."
    )
    parser.add_argument(
        "--region",
        type=str,
        default=AWS_REGION,
        help=f"AWS region for Bedrock. Default: {AWS_REGION}"
    )
    parser.add_argument(
        "--max-tokens",
        type=int,
        default=DEFAULT_MAX_TOKENS,
        help=f"Maximum tokens for the model's response. Default: {DEFAULT_MAX_TOKENS}"
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=DEFAULT_TEMPERATURE,
        help=f"Temperature for model's response generation. Default: {DEFAULT_TEMPERATURE}"
    )
    parser.add_argument(
        "--top-p",
        type=float,
        default=DEFAULT_TOP_P,
        help=f"Top P for model's response generation. Default: {DEFAULT_TOP_P}"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print full API response details."
    )

    args = parser.parse_args()

    try:
        # Initialize Bedrock client
        bedrock_client = boto3.client("bedrock-runtime", region_name=args.region)
        
        # Prepare the message payload
        messages = [
            {
                "role": "user",
                "content": [
                    {
                        "text": args.prompt
                    }
                ]
            }
        ]

        inference_config = {
            "maxTokens": args.max_tokens,
            "temperature": args.temperature,
            "topP": args.top_p
        }

        print(f"\nTesting Bedrock model: {args.model_id}")
        print(f"Prompt: \"{args.prompt}\"")
        print("Sending request to Amazon Bedrock...")
        
        # Send request to the model
        response = bedrock_client.converse(
            modelId=args.model_id,
            messages=messages,
            inferenceConfig=inference_config
        )

        # Extract and print the response
        response_text = response['output']['message']['content'][0]['text']
        print("\n[Model Response]")
        print(response_text)

        # Print full response if verbose mode is enabled
        if args.verbose:
            print("\n[Full API Response]")
            print(json.dumps(response, indent=2))

    except Exception as e:
        print(f"\nAn error occurred: {e}", file=sys.stderr)
        # Print more details from Bedrock error if available
        if hasattr(e, 'response') and 'Error' in e.response:
            print(f"Bedrock Error: {e.response['Error']}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
