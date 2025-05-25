#!/usr/bin/env python3

import argparse
import boto3
import os
import json
import uuid
import re
import sys
from botocore.exceptions import ClientError

# --- Configuration Constants ---
# Default Model ID for Amazon Titan
DEFAULT_BEDROCK_MODEL_ID = "us.amazon.nova-pro-v1:0"
# AWS Region for Bedrock and S3 services
AWS_REGION = "us-east-1" # You can change this if needed

# Document processing limits
MAX_BYTES_PAYLOAD_SIZE = 25 * 1024 * 1024  # 25MB in bytes
MAX_DOCS_FOR_BYTES_PAYLOAD = 5
MAX_DOCS_FOR_S3_PAYLOAD = 1000

# S3 configuration (if used)
S3_UPLOAD_PREFIX = "askpdf_temp_uploads/"

# Bedrock inference parameters
DEFAULT_INPUT_TOKENS = 30000  # Placeholder for typical model context window
DEFAULT_OUTPUT_TOKENS = 2048
DEFAULT_TOP_P = 0.9
DEFAULT_TEMPERATURE = 0.2


def sanitize_document_name(filename):
    """
    Sanitizes a filename to be a valid document name for Bedrock.
    Keeps alphanumeric, hyphens, parentheses, square brackets. Replaces others with underscore.
    Limits length to 60 characters.
    """
    base_name = os.path.splitext(filename)[0]
    sanitized = re.sub(r'[^a-zA-Z0-9\-\(\)\[\]_]', '_', base_name)
    return sanitized[:60]


def get_pdf_files_details(input_path):
    """
    Collects PDF file paths and their total size from a given file or directory.
    Returns a list of full file paths and the total size in bytes.
    """
    pdf_file_paths = []
    total_size_bytes = 0

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Error: Input path '{input_path}' not found.")

    if os.path.isfile(input_path):
        if input_path.lower().endswith(".pdf"):
            pdf_file_paths.append(input_path)
            total_size_bytes = os.path.getsize(input_path)
        else:
            raise ValueError(f"Error: Specified file '{input_path}' is not a PDF.")
    elif os.path.isdir(input_path):
        for root, _, files in os.walk(input_path):
            for file in files:
                if file.lower().endswith(".pdf"):
                    full_path = os.path.join(root, file)
                    pdf_file_paths.append(full_path)
                    total_size_bytes += os.path.getsize(full_path)
    else:
        raise ValueError(f"Error: Input path '{input_path}' is not a valid file or directory.")

    if not pdf_file_paths:
        raise FileNotFoundError(f"Error: No PDF files found at '{input_path}'.")

    return pdf_file_paths, total_size_bytes


def upload_to_s3(s3_client, bucket_name, file_path, s3_key):
    """Uploads a single file to S3."""
    try:
        s3_client.upload_file(file_path, bucket_name, s3_key)
        print(f"Successfully uploaded {os.path.basename(file_path)} to s3://{bucket_name}/{s3_key}")
    except Exception as e:
        raise RuntimeError(f"Error uploading {file_path} to S3 bucket {bucket_name}: {e}")


def delete_s3_objects(s3_client, bucket_name, s3_keys):
    """Deletes multiple objects from S3."""
    if not s3_keys:
        return
    objects_to_delete = {'Objects': [{'Key': key} for key in s3_keys]}
    try:
        s3_client.delete_objects(Bucket=bucket_name, Delete=objects_to_delete)
        print(f"Successfully deleted {len(s3_keys)} temporary objects from S3 bucket '{bucket_name}'.")
    except Exception as e:
        print(f"Warning: Failed to delete objects from S3 bucket '{bucket_name}': {e}. Manual cleanup may be required for keys: {s3_keys}")


def main():
    parser = argparse.ArgumentParser(
        description="Ask questions about PDF documents using Amazon Bedrock (Nova Pro model)."
    )
    parser.add_argument(
        "path",
        type=str,
        help="Path to a PDF file or a folder containing PDF files."
    )
    parser.add_argument(
        "question",
        type=str,
        help="The question to ask about the PDF(s)."
    )
    parser.add_argument(
        "--bucket",
        type=str,
        help="S3 bucket name for temporary storage if total PDF size > 25MB. Required in that case."
    )
    parser.add_argument(
        "--region",
        type=str,
        default=AWS_REGION,
        help=f"AWS region for Bedrock and S3. Default: {AWS_REGION}"
    )
    parser.add_argument(
        "--profile",
        type=str,
        help="AWS profile name to use for credentials"
    )
    parser.add_argument(
        "--output-tokens",
        type=int,
        default=DEFAULT_OUTPUT_TOKENS,
        help=f"Maximum tokens for the model's output response. Default: {DEFAULT_OUTPUT_TOKENS}"
    )
    parser.add_argument(
        "--input-tokens",
        type=int,
        default=DEFAULT_INPUT_TOKENS,
        help=f"Informational: Target/assumed input token capacity. Not directly enforced by the API for document inputs. Default: {DEFAULT_INPUT_TOKENS}"
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
        "--model-id",
        type=str,
        default=DEFAULT_BEDROCK_MODEL_ID,
        help=f"The Bedrock model ID to use. Default: {DEFAULT_BEDROCK_MODEL_ID}"
    )

    args = parser.parse_args()
    current_model_id = args.model_id

    try:
        pdf_file_paths, total_size_bytes = get_pdf_files_details(args.path)
    except (FileNotFoundError, ValueError) as e:
        print(e, file=sys.stderr)
        sys.exit(1)

    num_documents = len(pdf_file_paths)
    print(f"Found {num_documents} PDF document(s), total size: {total_size_bytes / (1024*1024):.2f} MB.")

    use_s3 = total_size_bytes > MAX_BYTES_PAYLOAD_SIZE
    s3_bucket_name = args.bucket
    s3_keys_uploaded = []

    # Initialize AWS session with optional profile
    session = boto3.Session(profile_name=args.profile) if args.profile else boto3.Session()
    bedrock_client = session.client("bedrock-runtime", region_name=args.region)
    s3_client = None
    aws_account_id = None

    if use_s3:
        if not s3_bucket_name:
            print(f"Error: Total PDF size ({total_size_bytes / (1024*1024):.2f} MB) exceeds {MAX_BYTES_PAYLOAD_SIZE / (1024*1024)} MB. "
                  f"Please provide an S3 bucket name using --bucket for temporary storage.", file=sys.stderr)
            sys.exit(1)
        if num_documents > MAX_DOCS_FOR_S3_PAYLOAD:
            print(f"Error: Number of documents ({num_documents}) exceeds the S3 limit of {MAX_DOCS_FOR_S3_PAYLOAD}.", file=sys.stderr)
            sys.exit(1)
        
        s3_client = session.client("s3", region_name=args.region)
        try:
            sts_client = session.client("sts", region_name=args.region)
            aws_account_id = sts_client.get_caller_identity()["Account"]
        except Exception as e:
            print(f"Error: Could not retrieve AWS Account ID for S3 bucket owner: {e}", file=sys.stderr)
            sys.exit(1)
        print(f"Using S3 bucket '{s3_bucket_name}' for document transfer.")
    else:
        if num_documents > MAX_DOCS_FOR_BYTES_PAYLOAD:
            print(f"Error: Number of documents ({num_documents}) exceeds the limit of {MAX_DOCS_FOR_BYTES_PAYLOAD} for direct upload. "
                  f"Total size is {total_size_bytes / (1024*1024):.2f} MB. "
                  f"If total size were > {MAX_BYTES_PAYLOAD_SIZE / (1024*1024)} MB, S3 would be used (up to {MAX_DOCS_FOR_S3_PAYLOAD} docs).", file=sys.stderr)
            sys.exit(1)
        print("Using direct byte transfer for documents.")


    document_content_list = []

    try:
        for idx, file_path in enumerate(pdf_file_paths):
            doc_name = sanitize_document_name(os.path.basename(file_path))
            unique_doc_name = f"doc_{idx}_{doc_name}" # Ensure unique and neutral

            if use_s3:
                s3_key = f"{S3_UPLOAD_PREFIX}{uuid.uuid4()}-{os.path.basename(file_path)}"
                upload_to_s3(s3_client, s3_bucket_name, file_path, s3_key)
                s3_keys_uploaded.append(s3_key)
                # For S3 documents, we need to download them first and then send as bytes
                # since the API doesn't directly support s3Location parameter
                s3_obj = s3_client.get_object(Bucket=s3_bucket_name, Key=s3_key)
                doc_bytes = s3_obj['Body'].read()
                document_source = {"bytes": doc_bytes}
            else: # Use bytes
                with open(file_path, "rb") as f:
                    doc_bytes = f.read()
                document_source = {"bytes": doc_bytes}

            document_content_list.append({
                "document": {
                    "format": "pdf",
                    "name": unique_doc_name,
                    "source": document_source
                }
            })

        # Add the question as the last part of the content
        content_payload = document_content_list + [{"text": args.question}]
        
        messages = [{"role": "user", "content": content_payload}]

        inference_config = {
            "maxTokens": args.output_tokens,
            "temperature": args.temperature,
            "topP": args.top_p
        }

        print(f"\nSending request to Amazon Bedrock ({current_model_id})...")
        response = bedrock_client.converse(
            modelId=current_model_id,
            messages=messages,
            inferenceConfig=inference_config
            # system_prompt could be added here if needed
        )

        response_text = response['output']['message']['content'][0]['text']
        print("\n[Model Response]")
        print(response_text)

        # Optional: Print full response for debugging
        # print("\n[Full API Response]")
        # print(json.dumps(response, indent=2))

    except ClientError as e:
        error_code = e.response.get("Error", {}).get("Code")
        error_message = e.response.get("Error", {}).get("Message", str(e))
        
        base_error_text = f"A Bedrock API error occurred with model {current_model_id}"
        specific_error_info = f"Bedrock Error Code: {error_code}, Message: {error_message}"

        full_message = f"\n{base_error_text}.\n{specific_error_info}"

        if error_code == "ValidationException" and "Input is too long" in error_message:
            full_message += (
                "\n\n[Additional Diagnostics]:"
                f"\nYour PDF file ({total_size_bytes / (1024*1024):.2f} MB) was uploaded successfully, "
                "but when Amazon Nova processed the document content (text, images, charts, etc.), "
                "it exceeded the model's token processing limit."
                f"\nModel used: {current_model_id}"
                "\n\nPossible solutions:"
                "\n1. Try a smaller document or split this PDF into smaller sections"
                "\n2. Try a different Nova model with a larger context window:"
                "\n   --model-id us.amazon.nova-premier-v1:0  (largest context window)"
                "\n   --model-id us.amazon.nova-lite-v1:0     (smaller but more efficient)"
                "\n3. The document may have many images/charts that consume significant tokens"
                "\n\nNote: File size (MB) ≠ token count. A PDF with lots of text/images can have "
                "high token usage even if the file size seems reasonable."
            )
        
        print(full_message, file=sys.stderr)
        sys.exit(1)
    except Exception as e: # General fallback for non-ClientError exceptions
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        if use_s3 and s3_keys_uploaded:
            print(f"\nCleaning up {len(s3_keys_uploaded)} S3 objects...")
            delete_s3_objects(s3_client, s3_bucket_name, s3_keys_uploaded)

if __name__ == "__main__":
    main()
