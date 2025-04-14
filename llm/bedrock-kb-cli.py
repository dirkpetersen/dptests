import argparse
import os
import boto3
import json
from botocore.exceptions import ClientError
import sys
from pathlib import Path

def create_role_if_needed():
    iam = boto3.client('iam')
    role_name = 'AmazonBedrockExecutionRoleForKnowledgeBase'
    try:
        return iam.get_role(RoleName=role_name)['Role']['Arn']
    except iam.exceptions.NoSuchEntityException:
        assume_role_policy = {
            "Version": "2012-10-17",
            "Statement": [{
                "Effect": "Allow",
                "Principal": {"Service": "bedrock.amazonaws.com"},
                "Action": "sts:AssumeRole"
            }]
        }
        role = iam.create_role(
            RoleName=role_name,
            AssumeRolePolicyDocument=json.dumps(assume_role_policy)
        )
        iam.attach_role_policy(
            RoleName=role_name,
            PolicyArn='arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess'
        )
        iam.attach_role_policy(
            RoleName=role_name,
            PolicyArn='arn:aws:iam::aws:policy/AWSBedrockFoundationModelPolicy'
        )
        return role['Role']['Arn']

def upload_to_s3(path, kb_name='default'):
    s3 = boto3.client('s3')
    bucket = 'my-bedrock-kb'
    prefix = kb_name + '/'
    
    # Clear existing prefix
    objects = s3.list_objects_v2(Bucket=bucket, Prefix=prefix).get('Contents', [])
    if objects:
        s3.delete_objects(Bucket=bucket, Delete={'Objects': [{'Key': obj['Key']} for obj in objects]})
    
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
    region = boto3.session.Session().region_name
    account_id = boto3.client('sts').get_caller_identity()['Account']
    
    try:
        bedrock.delete_knowledge_base(knowledgeBaseId=kb_name)
    except bedrock.exceptions.ResourceNotFoundException:
        pass
    
    role_arn = create_role_if_needed()
    
    return bedrock.create_knowledge_base(
        name=kb_name,
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
                "collectionArn": f"arn:aws:aoss:{region}:{account_id}:collection/your-collection-id",
                "fieldMapping": {
                    "vectorField": "bedrock-knowledge-base-default-vector",
                    "textField": "AMAZON_BEDROCK_TEXT_CHUNK",
                    "metadataField": "AMAZON_BEDROCK_METADATA"
                }
            }
        }
    )

def ask_question(kb_name='default', query=None):
    bedrock = boto3.client('bedrock-agent-runtime')
    while True:
        try:
            response = bedrock.retrieve_and_generate(
                input={'text': query or input("\nQuestion: ")},
                retrieveAndGenerateConfiguration={
                    "type": "KNOWLEDGE_BASE",
                    "knowledgeBaseConfiguration": {
                        "knowledgeBaseId": kb_name,
                        "modelArn": "arn:aws:bedrock:us-west-2::foundation-model/anthropic.claude-3-sonnet-20240229-v1:0"
                    }
                }
            )
            print(f"\n**Answer**:\n{response['output']['text']}\n")
            if query: break
        except KeyboardInterrupt:
            print("\nExiting...")
            sys.exit(0)

def main():
    parser = argparse.ArgumentParser(description='Bedrock Knowledge Base CLI')
    subparsers = parser.add_subparsers(dest='command')
    
    # Upload command
    upload_parser = subparsers.add_parser('upload')
    upload_parser.add_argument('path', help='File or directory path to upload')
    upload_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    # Ask command
    ask_parser = subparsers.add_parser('ask')
    ask_parser.add_argument('question', nargs='?', help='Question to ask (optional)')
    ask_parser.add_argument('--kb', default='default', help='Knowledge base name')
    
    args = parser.parse_args()
    
    if args.command == 'upload':
        upload_to_s3(args.path, args.kb)
        create_knowledge_base(args.kb)
        print(f"Knowledge base '{args.kb}' updated!")
    elif args.command == 'ask':
        ask_question(args.kb, args.question)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
