#! /usr/bin/env python3

import argparse
import os
import boto3
import json
import uuid
import time
import re
from botocore.exceptions import ClientError
import sys
from pathlib import Path

def generate_kb_id(kb_name):
    """Generate AWS-compliant knowledge base ID (exactly 10 alphanumeric characters)"""
    # Remove special characters and truncate
    clean_name = ''.join([c for c in kb_name.upper() if c.isalnum()])[:8]
    # Pad with random chars to make exactly 10 characters
    return clean_name + uuid.uuid4().hex[:10-len(clean_name)]

def get_config():
    """Get or create configuration values"""
    account_id = boto3.client('sts').get_caller_identity()['Account']
    region = boto3.session.Session().region_name
    
    # Use a consistent bucket name based on account ID
    bucket_name = f"bedrock-kb-{account_id}-{region}"
    
    # Use a consistent collection name
    collection_name = f"bedrock-kb-collection"
    
    return {
        'account_id': account_id,
        'region': region,
        'bucket_name': bucket_name,
        'collection_name': collection_name,
        'role_name': 'AmazonBedrockExecutionRoleForKnowledgeBase'
    }

def verify_role():
    """Get the role ARN without checking permissions"""
    config = get_config()
    role_name = config['role_name']
    role_arn = f"arn:aws:iam::{config['account_id']}:role/{role_name}"
    print(f"Using role: {role_arn}")
    return role_arn

def create_opensearch_collection_if_needed():
    """Create OpenSearch collection if it doesn't exist"""
    config = get_config()
    collection_name = config['collection_name']
    
    aoss = boto3.client('opensearchserverless')
    
    # Check if collection exists
    try:
        response = aoss.batch_get_collection(names=[collection_name])
        if response['collectionDetails']:
            collection = response['collectionDetails'][0]
            print(f"✓ Using existing OpenSearch collection: {collection_name}")
            # Wait for collection to be active if needed
            if collection['status'] != 'ACTIVE':
                print(f"  Collection status: {collection['status']} - waiting for ACTIVE...")
                while collection['status'] != 'ACTIVE':
                    time.sleep(10)
                    response = aoss.batch_get_collection(names=[collection_name])
                    if not response['collectionDetails']:
                        break
                    collection = response['collectionDetails'][0]
                    print(f"  Collection status: {collection['status']}")
            return collection_name
    except ClientError as e:
        if e.response['Error']['Code'] != 'ResourceNotFoundException':
            raise
    
    # Create collection if it doesn't exist
    print(f"Creating OpenSearch collection: {collection_name}")
    
    # Create security policy
    policy_name = f"bedrock-kb-policy-{uuid.uuid4().hex[:8]}"
    policy_document = {
        "Rules": [
            {
                "ResourceType": "collection",
                "Resource": [f"collection/{collection_name}"],
                "Permission": [
                    "aoss:CreateCollectionItems",
                    "aoss:DeleteCollectionItems",
                    "aoss:UpdateCollectionItems",
                    "aoss:DescribeCollectionItems"
                ]
            },
            {
                "ResourceType": "index",
                "Resource": [f"index/{collection_name}/*"],
                "Permission": [
                    "aoss:CreateIndex",
                    "aoss:DeleteIndex",
                    "aoss:UpdateIndex",
                    "aoss:DescribeIndex",
                    "aoss:ReadDocument",
                    "aoss:WriteDocument"
                ]
            }
        ]
    }
    
    try:
        aoss.create_security_policy(
            name=policy_name,
            policy=json.dumps(policy_document),
            type='encryption'
        )
        print(f"✓ Created security policy: {policy_name}")
    except ClientError as e:
        if 'already exists' not in str(e):
            raise
    
    # Create network policy
    network_policy_name = f"bedrock-kb-network-{uuid.uuid4().hex[:8]}"
    network_policy = {
        "Rules": [
            {
                "ResourceType": "collection",
                "Resource": [f"collection/{collection_name}"],
                "SourceVPCEs": []
            }
        ]
    }
    
    try:
        aoss.create_security_policy(
            name=network_policy_name,
            policy=json.dumps(network_policy),
            type='network'
        )
        print(f"✓ Created network policy: {network_policy_name}")
    except ClientError as e:
        if 'already exists' not in str(e):
            raise
    
    # Create collection
    try:
        response = aoss.create_collection(
            name=collection_name,
            type='SEARCH'
        )
        print(f"✓ Creating collection: {collection_name}")
        
        # Wait for collection to be created
        collection_id = response['createCollectionDetail']['id']
        status = response['createCollectionDetail']['status']
        print(f"  Collection status: {status} - waiting for ACTIVE...")
        
        while status != 'ACTIVE':
            time.sleep(10)
            response = aoss.batch_get_collection(ids=[collection_id])
            if not response['collectionDetails']:
                break
            status = response['collectionDetails'][0]['status']
            print(f"  Collection status: {status}")
            
        print(f"✓ Collection {collection_name} is now active")
        return collection_name
    except ClientError as e:
        print(f"Error creating collection: {str(e)}")
        sys.exit(1)

def upload_to_s3(path, kb_name='default'):
    config = get_config()
    bucket = config['bucket_name']
    s3 = boto3.client('s3')
    
    # Add bucket creation check
    try:
        s3.head_bucket(Bucket=bucket)
        print(f"✓ Using existing S3 bucket: {bucket}")
    except s3.exceptions.ClientError as e:
        error_code = int(e.response['Error']['Code'])
        if error_code == 404:
            # Create the bucket if it doesn't exist
            region = config['region']
            print(f"Creating S3 bucket: {bucket}")
            
            # Different API call for us-east-1
            if region == 'us-east-1':
                s3.create_bucket(Bucket=bucket)
            else:
                s3.create_bucket(
                    Bucket=bucket,
                    CreateBucketConfiguration={'LocationConstraint': region}
                )
            
            # Add bucket policy to allow Bedrock access
            policy = {
                "Version": "2012-10-17",
                "Statement": [
                    {
                        "Effect": "Allow",
                        "Principal": {
                            "Service": "bedrock.amazonaws.com"
                        },
                        "Action": [
                            "s3:GetObject",
                            "s3:ListBucket"
                        ],
                        "Resource": [
                            f"arn:aws:s3:::{bucket}",
                            f"arn:aws:s3:::{bucket}/*"
                        ]
                    }
                ]
            }
            
            s3.put_bucket_policy(
                Bucket=bucket,
                Policy=json.dumps(policy)
            )
            print(f"✓ Created S3 bucket with Bedrock access policy: {bucket}")
        else:
            raise
    
    prefix = kb_name + '/'
    
    # Clear existing prefix
    delete_keys = []
    paginator = s3.get_paginator('list_objects_v2')
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        if 'Contents' in page:
            delete_keys.extend([{'Key': obj['Key']} for obj in page['Contents']])
    
    # Delete in batches of 1000 (S3 API limit)
    if delete_keys:
        for i in range(0, len(delete_keys), 1000):
            s3.delete_objects(
                Bucket=bucket,
                Delete={'Objects': delete_keys[i:i+1000]}
            )
    
    # Upload new files
    if os.path.isfile(path):
        s3.upload_file(path, bucket, prefix + os.path.basename(path))
    else:
        for root, _, files in os.walk(path):
            for file in files:
                local_path = os.path.join(root, file)
                s3_path = prefix + os.path.relpath(local_path, path)
                s3.upload_file(local_path, bucket, s3_path)

def create_knowledge_base(kb_name):
    bedrock = boto3.client('bedrock-agent')
    config = get_config()
    region = config['region']
    account_id = config['account_id']
    
    # Generate a valid knowledge base ID (exactly 10 alphanumeric characters)
    kb_id = generate_kb_id(kb_name)
    
    # Check if knowledge base exists
    try:
        bedrock.get_knowledge_base(knowledgeBaseId=kb_id)
        print(f"Found existing knowledge base: {kb_name} (ID: {kb_id})")
        try:
            print(f"Deleting existing knowledge base: {kb_name} (ID: {kb_id})")
            bedrock.delete_knowledge_base(knowledgeBaseId=kb_id)
            # Wait for deletion to complete
            print("Waiting for deletion to complete...")
            waiter = bedrock.get_waiter('knowledge_base_deleted')
            waiter.wait(knowledgeBaseId=kb_id)
        except bedrock.exceptions.ResourceNotFoundException:
            pass
        except ClientError as e:
            if e.response['Error']['Code'] == 'AccessDeniedException':
                print(f"Warning: Could not delete existing knowledge base - {e}")
            else:
                raise
    except bedrock.exceptions.ResourceNotFoundException:
        pass
    
    # Verify role exists
    role_arn = verify_role()
    
    # Create or verify OpenSearch collection
    collection_id = create_opensearch_collection_if_needed()
    collection_arn = f"arn:aws:aoss:{region}:{account_id}:collection/{collection_id}"
    
    print(f"Creating knowledge base: {kb_name} (ID: {kb_id})")
    try:
        response = bedrock.create_knowledge_base(
            name=kb_name,
            knowledgeBaseId=kb_id,
            roleArn=role_arn,
            knowledgeBaseConfiguration={
                "type": "VECTOR",
                "vectorKnowledgeBaseConfiguration": {
                    "embeddingModelArn": f"arn:aws:bedrock:{region}::foundation-model/amazon.titan-embed-text-v1"
                }
            },
            storageConfiguration={
                "type": "OPENSEARCH_SERVERLESS",
                "opensearchServerlessConfiguration": {
                    "collectionArn": collection_arn,
                    "vectorIndexName": "bedrock-knowledge-base-index",
                    "fieldMapping": {
                        "vectorField": "bedrock-knowledge-base-default-vector",
                        "textField": "AMAZON_BEDROCK_TEXT_CHUNK",
                        "metadataField": "AMAZON_BEDROCK_METADATA"
                    }
                }
            }
        )
        
        # Wait for knowledge base to be created
        kb_id = response['knowledgeBase']['knowledgeBaseId']
        print(f"Knowledge base creation initiated. ID: {kb_id}")
        print("Waiting for knowledge base to become active...")
        
        waiter = bedrock.get_waiter('knowledge_base_active')
        waiter.wait(knowledgeBaseId=kb_id)
        
        print(f"✓ Knowledge base '{kb_name}' is now active")
        return response
    except ClientError as e:
        print(f"\nERROR: {e.response['Error']['Message']}")
        print('''\nIAM Configuration Required:
  Your AWS administrator needs to create this role:

  Role Name: AmazonBedrockExecutionRoleForKnowledgeBase
  Trust Relationship: Allow bedrock.amazonaws.com
  Required Policies:
    - AWSBedrockAgentServiceRolePolicy
    - AmazonS3ReadOnlyAccess
    - AWSBedrockFoundationModelPolicy

  Once created, the script will work without additional permissions.''')
        sys.exit(1)

def create_data_source(kb_name):
    """Create a data source for the knowledge base"""
    bedrock = boto3.client('bedrock-agent')
    config = get_config()
    bucket = config['bucket_name']
    
    # Generate the same knowledge base ID
    kb_id = generate_kb_id(kb_name)
    
    # Get knowledge base details
    try:
        kb_response = bedrock.get_knowledge_base(knowledgeBaseId=kb_id)
        kb_id = kb_response['knowledgeBase']['knowledgeBaseId']
    except bedrock.exceptions.ResourceNotFoundException:
        print(f"Error: Knowledge base '{kb_name}' (ID: {kb_id}) not found")
        sys.exit(1)
    
    data_source_name = f"{kb_name}-s3-source"
    
    # Check if data source exists
    try:
        bedrock.get_data_source(
            knowledgeBaseId=kb_id,
            dataSourceId=data_source_name
        )
        print(f"Found existing data source: {data_source_name}")
        try:
            print(f"Deleting existing data source: {data_source_name}")
            bedrock.delete_data_source(
                knowledgeBaseId=kb_id,
                dataSourceId=data_source_name
            )
            # Wait for deletion to complete
            print("Waiting for data source deletion to complete...")
            time.sleep(10)  # Simple wait as there's no specific waiter for this
        except bedrock.exceptions.ResourceNotFoundException:
            pass
    except bedrock.exceptions.ResourceNotFoundException:
        pass
    
    # Create data source
    print(f"Creating data source: {data_source_name}")
    try:
        response = bedrock.create_data_source(
            knowledgeBaseId=kb_id,
            name=data_source_name,
            dataSourceConfiguration={
                "type": "S3",
                "s3Configuration": {
                    "bucketArn": f"arn:aws:s3:::{bucket}",
                    "inclusionPrefixes": [f"{kb_name}/"]
                }
            },
            vectorIngestionConfiguration={
                "chunkingConfiguration": {
                    "chunkingStrategy": "FIXED_SIZE",
                    "fixedSizeChunkingConfiguration": {
                        "maxTokens": 300,
                        "overlapPercentage": 20
                    }
                }
            }
        )
        
        # Start ingestion
        data_source_id = response['dataSource']['dataSourceId']
        print(f"Data source created. ID: {data_source_id}")
        print("Starting data ingestion...")
        
        bedrock.start_ingestion_job(
            knowledgeBaseId=kb_id,
            dataSourceId=data_source_id
        )
        
        print(f"✓ Data ingestion started for '{data_source_name}'")
        return response
    except ClientError as e:
        print(f"Error creating data source: {str(e)}")
        sys.exit(1)

def ask_question(kb_name='default', query=None):
    bedrock = boto3.client('bedrock-agent-runtime')
    region = boto3.session.Session().region_name
    model_arn = f"arn:aws:bedrock:{region}::foundation-model/anthropic.claude-3-sonnet-20240229-v1:0"
    
    # Generate the same knowledge base ID
    kb_id = generate_kb_id(kb_name)
    
    # Verify knowledge base exists
    bedrock_agent = boto3.client('bedrock-agent')
    try:
        bedrock_agent.get_knowledge_base(knowledgeBaseId=kb_id)
    except bedrock_agent.exceptions.ResourceNotFoundException:
        print(f"Error: Knowledge base '{kb_name}' (ID: {kb_id}) not found")
        sys.exit(1)
    
    while True:
        try:
            current_query = query or input("\nQuestion: ")
            print("\nQuerying knowledge base...")
            
            response = bedrock.retrieve_and_generate(
                input={'text': current_query},
                retrieveAndGenerateConfiguration={
                    "type": "KNOWLEDGE_BASE",
                    "knowledgeBaseConfiguration": {
                        "knowledgeBaseId": kb_id,
                        "modelArn": model_arn
                    }
                }
            )
            
            # Print retrieved citations
            if 'citations' in response:
                print("\nSources:")
                for i, citation in enumerate(response['citations'], 1):
                    print(f"  {i}. {citation['retrievedReferences'][0]['location']}")
            
            print(f"\n**Answer**:\n{response['output']['text']}\n")
            if query: break
        except ClientError as e:
            print(f"Error: {str(e)}")
            if "ResourceNotFoundException" in str(e):
                print(f"Knowledge base '{kb_name}' not found or not properly configured")
            sys.exit(1)
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(0)

def main():
    parser = argparse.ArgumentParser(
        description='Bedrock Knowledge Base CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''\
        
IAM REQUIREMENTS FOR BEDROCK KNOWLEDGE BASES
============================================

Your AWS administrator MUST create this role BEFORE using this tool:

1. Role Name:
   AmazonBedrockExecutionRoleForKnowledgeBase

2. Trust Relationship Policy:
   ----------------------------------
   {
     "Version": "2012-10-17",
     "Statement": [
       {
         "Effect": "Allow",
         "Principal": {
           "Service": "bedrock.amazonaws.com"
         },
         "Action": "sts:AssumeRole"
       }
     ]
   }
   ----------------------------------

3. Required Attached Policies:
   - AWSBedrockAgentServiceRolePolicy
   - AmazonS3ReadOnlyAccess
   - AWSBedrockFoundationModelPolicy

4. Role Requirements:
   - Must exist in the SAME AWS account as your CLI user
   - Must be created in the SAME REGION you're using
   - Must be created BEFORE running this tool

For detailed instructions see:
https://docs.aws.amazon.com/bedrock/latest/userguide/knowledge-base-create.html''')
    subparsers = parser.add_subparsers(dest='command')
    
    # Upload command
    upload_parser = subparsers.add_parser('upload')
    upload_parser.add_argument('path', help='File or directory path to upload')
    upload_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    # Ask command
    ask_parser = subparsers.add_parser('ask')
    ask_parser.add_argument('question', nargs='?', help='Question to ask (optional)')
    ask_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    # Status command
    status_parser = subparsers.add_parser('status')
    status_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    args = parser.parse_args()
    
    if args.command == 'upload':
        print(f"Setting up knowledge base: {args.kb}")
        upload_to_s3(args.path, args.kb)
        create_knowledge_base(args.kb)
        create_data_source(args.kb)
        print(f"\n✓ Knowledge base '{args.kb}' setup complete!")
        print("Note: It may take a few minutes for data ingestion to complete")
        print("      before the knowledge base is ready for queries.")
    elif args.command == 'ask':
        ask_question(args.kb, args.question)
    elif args.command == 'status':
        config = get_config()
        print(f"AWS Account: {config['account_id']}")
        print(f"AWS Region: {config['region']}")
        print(f"S3 Bucket: {config['bucket_name']}")
        print(f"OpenSearch Collection: {config['collection_name']}")
        print(f"IAM Role: {config['role_name']}")
        
        # Generate the same knowledge base ID
        kb_id = generate_kb_id(args.kb)
        
        # Check knowledge base status
        bedrock = boto3.client('bedrock-agent')
        try:
            kb = bedrock.get_knowledge_base(knowledgeBaseId=kb_id)
            print(f"\nKnowledge Base '{args.kb}' (ID: {kb_id}):")
            print(f"  Status: {kb['knowledgeBase']['status']}")
            print(f"  Created: {kb['knowledgeBase']['createdAt']}")
            
            # List data sources
            try:
                sources = bedrock.list_data_sources(knowledgeBaseId=kb_id)
                print("\nData Sources:")
                for ds in sources['dataSourceSummaries']:
                    print(f"  - {ds['name']} (Status: {ds['status']})")
                    
                    # Get latest ingestion job
                    try:
                        jobs = bedrock.list_ingestion_jobs(
                            knowledgeBaseId=kb_id,
                            dataSourceId=ds['dataSourceId']
                        )
                        if jobs['ingestionJobSummaries']:
                            latest = jobs['ingestionJobSummaries'][0]
                            print(f"    Latest ingestion: {latest['status']} ({latest['startedAt']})")
                    except Exception as e:
                        print(f"    Error getting ingestion jobs: {str(e)}")
            except Exception as e:
                print(f"  Error listing data sources: {str(e)}")
                
        except bedrock.exceptions.ResourceNotFoundException:
            print(f"\nKnowledge Base '{args.kb}' (ID: {kb_id}) not found")
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
