#! /usr/bin/env python3

#python3 -m pip install huggingface_hub boto3

from huggingface_hub import snapshot_download
import boto3
import os

model = "DeepSeek-R1-Distill-Llama-70B"


# Example using the 8B distilled model
model_id = f"deepseek-ai/{model}"
local_dir = snapshot_download(
    repo_id = model_id, 
    local_dir = model
)

s3_client = boto3.client('s3', region_name='us-west-2')
bucket_name = 'arcs-temp'
local_directory = model

# Upload all model files to S3
for root, dirs, files in os.walk(local_directory):
    for file in files:
        local_path = os.path.join(root, file)
        s3_key = os.path.relpath(local_path, local_directory)
        print(f'{model}/{s3_key}')
        s3_client.upload_file(local_path, bucket_name, f'{model}/{s3_key}')


