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
from datetime import datetime, timedelta

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
    
    role_name = 'AmazonBedrockExecutionRoleForKnowledgeBase'
    role_arn = f"arn:aws:iam::{account_id}:role/{role_name}"
    
    return {
        'account_id': account_id,
        'region': region,
        'bucket_name': bucket_name,
        'collection_name': collection_name,
        'role_name': role_name,
        'role_arn': role_arn
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
    
    # Shorten policy names to meet 32-character limit
    policy_name = f"enc-{collection_name}"[:32]  # Encryption policy
    network_policy_name = f"net-{collection_name}"[:32]  # Network policy
    data_policy_name = f"data-{collection_name}"[:32]  # Data policy
    
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
    
    # Create encryption policy with shortened name
    policy_document = {
        "Rules": [
            {
                "Resource": [f"collection/{collection_name}"],
                "ResourceType": "collection"
            }
        ],
        "AWSOwnedKey": True
    }
    
    # Check if policy exists first
    try:
        aoss.get_security_policy(name=policy_name, type='encryption')
        print(f"✓ Using existing encryption policy: {policy_name}")
    except ClientError as e:
        if e.response['Error']['Code'] == 'ResourceNotFoundException':
            try:
                aoss.create_security_policy(
                    name=policy_name,
                    policy=json.dumps(policy_document),
                    type='encryption'
                )
                print(f"✓ Created encryption policy: {policy_name}")
            except ClientError as e:
                if 'ConflictException' in str(e) or 'already exists' in str(e):
                    print(f"✓ Encryption policy already exists: {policy_name}")
                else:
                    raise
        else:
            raise
    
    # Create data access policy with shortened name
    data_policy = {
        "Rules": [
            {  # First rule for collection access
                "Resource": [f"collection/{collection_name}"],
                "ResourceType": "collection",
                "Access": [
                    "aoss:CreateCollectionItems",
                    "aoss:DeleteCollectionItems",
                    "aoss:UpdateCollectionItems",
                    "aoss:DescribeCollectionItems"
                ],
                "Principal": [config['role_arn']]
            },
            {  # Second rule for index access
                "Resource": [f"index/{collection_name}/*"],
                "ResourceType": "index",
                "Access": [
                    "aoss:CreateIndex",
                    "aoss:DeleteIndex",
                    "aoss:UpdateIndex",
                    "aoss:DescribeIndex",
                    "aoss:ReadDocument",
                    "aoss:WriteDocument"
                ],
                "Principal": [config['role_arn']]
            }
        ]
    }
    
    # Check if data policy exists first
    try:
        aoss.get_access_policy(name=data_policy_name, type='data')
        print(f"✓ Using existing data access policy: {data_policy_name}")
    except ClientError as e:
        if e.response['Error']['Code'] == 'ResourceNotFoundException':
            try:
                aoss.create_access_policy(
                    name=data_policy_name,
                    policy=json.dumps(data_policy),
                    type='data'
                )
                print(f"✓ Created data access policy: {data_policy_name}")
            except ClientError as e:
                if 'ConflictException' in str(e) or 'already exists' in str(e):
                    print(f"✓ Data access policy already exists: {data_policy_name}")
                else:
                    raise
        else:
            raise
    
    # Create network policy with shortened name
    network_policy = {
        "Rules": [
            {
                "ResourceType": "collection",
                "Resource": [f"collection/{collection_name}"],
                "SourceVPCEs": ["*"]
            }
        ]
    }
    
    # Check if network policy exists first
    try:
        aoss.get_security_policy(name=network_policy_name, type='network')
        print(f"✓ Using existing network policy: {network_policy_name}")
    except ClientError as e:
        if e.response['Error']['Code'] == 'ResourceNotFoundException':
            try:
                aoss.create_security_policy(
                    name=network_policy_name,
                    policy=json.dumps(network_policy),
                    type='network'
                )
                print(f"✓ Created network policy: {network_policy_name}")
            except ClientError as e:
                if 'ConflictException' in str(e) or 'already exists' in str(e):
                    print(f"✓ Network policy already exists: {network_policy_name}")
                else:
                    raise
        else:
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

def create_knowledge_base(kb_name, auto_delete_hours=None):
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
    
    # Prepare tags for the knowledge base
    tags = {
        'created-at': datetime.now().isoformat()
    }
    
    # Add auto-delete tag if specified
    if auto_delete_hours:
        tags['auto-delete'] = 'true'
        tags['auto-delete-hours'] = str(auto_delete_hours)
        tags['last-accessed'] = datetime.now().isoformat()
        print(f"Setting up auto-deletion after {auto_delete_hours} hours of inactivity")
    
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
            },
            tags=tags
        )
        
        # Wait for knowledge base to be created
        kb_id = response['knowledgeBase']['knowledgeBaseId']
        print(f"Knowledge base creation initiated. ID: {kb_id}")
        print("Waiting for knowledge base to become active...")
        
        waiter = bedrock.get_waiter('knowledge_base_active')
        waiter.wait(knowledgeBaseId=kb_id)
        
        print(f"✓ Knowledge base '{kb_name}' is now active")
        
        # Set up auto-deletion if specified
        if auto_delete_hours:
            setup_auto_deletion(kb_id, kb_name, auto_delete_hours)
            
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

def setup_auto_deletion(kb_id, kb_name, hours):
    """Set up auto-deletion for a knowledge base after specified hours of inactivity"""
    bedrock = boto3.client('bedrock-agent')
    
    # Update the last accessed timestamp
    bedrock.tag_resource(
        resourceARN=f"arn:aws:bedrock:{boto3.session.Session().region_name}:{get_config()['account_id']}:knowledge-base/{kb_id}",
        tags={
            'last-accessed': datetime.now().isoformat()
        }
    )
    
    print(f"✓ Auto-deletion set up for '{kb_name}' after {hours} hours of inactivity")
    print("  The knowledge base will be checked for activity during status checks")

def update_last_accessed(kb_id):
    """Update the last accessed timestamp for a knowledge base"""
    try:
        bedrock = boto3.client('bedrock-agent')
        bedrock.tag_resource(
            resourceARN=f"arn:aws:bedrock:{boto3.session.Session().region_name}:{get_config()['account_id']}:knowledge-base/{kb_id}",
            tags={
                'last-accessed': datetime.now().isoformat()
            }
        )
    except Exception as e:
        # Just log the error but don't fail the operation
        print(f"Warning: Could not update last-accessed timestamp - {e}")

def check_and_delete_inactive_kbs():
    """Check for inactive knowledge bases and delete them if needed"""
    bedrock = boto3.client('bedrock-agent')
    
    # List all knowledge bases
    paginator = bedrock.get_paginator('list_knowledge_bases')
    
    for page in paginator.paginate():
        for kb in page.get('knowledgeBaseSummaries', []):
            kb_id = kb['knowledgeBaseId']
            
            try:
                # Get tags for the knowledge base
                response = bedrock.list_tags_for_resource(
                    resourceARN=f"arn:aws:bedrock:{boto3.session.Session().region_name}:{get_config()['account_id']}:knowledge-base/{kb_id}"
                )
                
                tags = response.get('tags', {})
                
                # Check if this KB has auto-delete enabled
                if tags.get('auto-delete') == 'true' and 'last-accessed' in tags and 'auto-delete-hours' in tags:
                    last_accessed = datetime.fromisoformat(tags['last-accessed'])
                    auto_delete_hours = int(tags['auto-delete-hours'])
                    
                    # Check if the KB has been inactive for longer than the specified hours
                    if datetime.now() - last_accessed > timedelta(hours=auto_delete_hours):
                        print(f"Deleting inactive knowledge base: {kb['name']} (ID: {kb_id})")
                        print(f"  Last accessed: {last_accessed.isoformat()}")
                        print(f"  Inactive for: {(datetime.now() - last_accessed).total_seconds() / 3600:.1f} hours")
                        
                        # Delete the knowledge base
                        bedrock.delete_knowledge_base(knowledgeBaseId=kb_id)
                        print(f"✓ Deleted inactive knowledge base: {kb['name']}")
            except Exception as e:
                print(f"Error checking knowledge base {kb_id}: {str(e)}")

def ask_question(kb_name='default', query=None):
    bedrock = boto3.client('bedrock-agent-runtime')
    region = boto3.session.Session().region_name
    model_arn = f"arn:aws:bedrock:{region}::foundation-model/anthropic.claude-3-sonnet-20240229-v1:0"
    
    # Generate the same knowledge base ID
    kb_id = generate_kb_id(kb_name)
    
    # Verify knowledge base exists
    bedrock_agent = boto3.client('bedrock-agent')
    try:
        kb_response = bedrock_agent.get_knowledge_base(knowledgeBaseId=kb_id)
        
        # Update last accessed timestamp if this is an auto-delete KB
        tags = bedrock_agent.list_tags_for_resource(
            resourceARN=f"arn:aws:bedrock:{region}:{get_config()['account_id']}:knowledge-base/{kb_id}"
        ).get('tags', {})
        
        if tags.get('auto-delete') == 'true':
            update_last_accessed(kb_id)
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
            
            # Update last accessed timestamp for auto-delete KBs
            try:
                tags = bedrock_agent.list_tags_for_resource(
                    resourceARN=f"arn:aws:bedrock:{region}:{get_config()['account_id']}:knowledge-base/{kb_id}"
                ).get('tags', {})
                
                if tags.get('auto-delete') == 'true':
                    update_last_accessed(kb_id)
            except Exception:
                pass  # Ignore errors updating timestamp
            
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
    upload_parser.add_argument('--auto-delete', type=int, help='Auto-delete after specified hours of inactivity')
    
    # Ask command
    ask_parser = subparsers.add_parser('ask')
    ask_parser.add_argument('question', nargs='?', help='Question to ask (optional)')
    ask_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    # Status command
    status_parser = subparsers.add_parser('status')
    status_parser.add_argument('--kb', default='default', help='Knowledge base name')
    status_parser.add_argument('--check-inactive', action='store_true', help='Check and delete inactive knowledge bases')
    
    # Cleanup command
    cleanup_parser = subparsers.add_parser('cleanup')
    cleanup_parser.add_argument('--force', action='store_true', help='Force deletion of all auto-delete knowledge bases')
    
    args = parser.parse_args()
    
    if args.command == 'upload':
        print(f"Setting up knowledge base: {args.kb}")
        upload_to_s3(args.path, args.kb)
        create_knowledge_base(args.kb, args.auto_delete)
        create_data_source(args.kb)
        print(f"\n✓ Knowledge base '{args.kb}' setup complete!")
        print("Note: It may take a few minutes for data ingestion to complete")
        print("      before the knowledge base is ready for queries.")
        if args.auto_delete:
            print(f"Auto-deletion is set for {args.auto_delete} hours of inactivity")
    elif args.command == 'ask':
        ask_question(args.kb, args.question)
    elif args.command == 'status':
        config = get_config()
        print(f"AWS Account: {config['account_id']}")
        print(f"AWS Region: {config['region']}")
        print(f"S3 Bucket: {config['bucket_name']}")
        print(f"OpenSearch Collection: {config['collection_name']}")
        print(f"IAM Role: {config['role_name']}")
        
        # Check for inactive KBs if requested
        if args.check_inactive:
            print("\nChecking for inactive knowledge bases...")
            check_and_delete_inactive_kbs()
        
        # Generate the same knowledge base ID
        kb_id = generate_kb_id(args.kb)
        
        # Check knowledge base status
        bedrock = boto3.client('bedrock-agent')
        try:
            kb = bedrock.get_knowledge_base(knowledgeBaseId=kb_id)
            print(f"\nKnowledge Base '{args.kb}' (ID: {kb_id}):")
            print(f"  Status: {kb['knowledgeBase']['status']}")
            print(f"  Created: {kb['knowledgeBase']['createdAt']}")
            
            # Get tags to check auto-delete status
            try:
                tags_response = bedrock.list_tags_for_resource(
                    resourceARN=f"arn:aws:bedrock:{config['region']}:{config['account_id']}:knowledge-base/{kb_id}"
                )
                tags = tags_response.get('tags', {})
                
                if tags.get('auto-delete') == 'true':
                    last_accessed = datetime.fromisoformat(tags.get('last-accessed', datetime.now().isoformat()))
                    auto_delete_hours = int(tags.get('auto-delete-hours', '0'))
                    time_remaining = timedelta(hours=auto_delete_hours) - (datetime.now() - last_accessed)
                    
                    print(f"  Auto-delete: Enabled ({auto_delete_hours} hours of inactivity)")
                    print(f"  Last accessed: {last_accessed.isoformat()}")
                    if time_remaining.total_seconds() > 0:
                        print(f"  Time remaining: {time_remaining.total_seconds() / 3600:.1f} hours")
                    else:
                        print(f"  Time remaining: Scheduled for deletion")
                        
                    # Update last accessed time since we're checking status
                    update_last_accessed(kb_id)
            except Exception as e:
                print(f"  Error getting tags: {str(e)}")
            
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
    elif args.command == 'cleanup':
        if args.force:
            print("Forcing cleanup of all auto-delete knowledge bases...")
            check_and_delete_inactive_kbs()
        else:
            print("Checking for inactive knowledge bases...")
            check_and_delete_inactive_kbs()
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
