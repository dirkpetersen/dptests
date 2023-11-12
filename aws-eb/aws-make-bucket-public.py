#! /usr/bin/env python3

import boto3
import json

# Initialize a boto3 client
s3_client = boto3.client('s3')

# Define the bucket name
bucket_name = 'easybuild-cache'

# Bucket policy that grants public read access
public_policy = {
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Principal": "*",
            "Action": ["s3:GetObject"],
            "Resource": [f"arn:aws:s3:::{bucket_name}/*"]
        }
    ]
}

# Update the bucket policy
s3_client.put_bucket_policy(
    Bucket=bucket_name,
    Policy=json.dumps(public_policy)
)

s3_client.put_bucket_request_payment(
    Bucket=bucket_name,
    RequestPaymentConfiguration={
        'Payer': 'Requester'
    }
)

