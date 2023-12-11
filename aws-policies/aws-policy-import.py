#! /usr/bin/env python3

import boto3, json, os, sys, argparse
from botocore.exceptions import ClientError

def create_policy(policy_name, policy_document):
    try:
        iam.create_policy(
            PolicyName=policy_name,
            PolicyDocument=json.dumps(policy_document)
        )
        print(f"Policy {policy_name} created successfully.")
    except ClientError as error:
        print(f"Failed to create policy {policy_name}: {error}")

def attach_policy_to_user(policy_arn, user_name):
    try:
        iam.attach_user_policy(
            UserName=user_name,
            PolicyArn=policy_arn
        )
        print(f"Attached policy {policy_arn} to user {user_name}.")
    except ClientError as error:
        print(f"Failed to attach policy {policy_arn} to user {user_name}: {error}")

# Parse command line argument
parser = argparse.ArgumentParser(description="Create IAM policies from JSON files in a folder.")
parser.add_argument("folder", help="Path to the folder containing JSON policy files")
args = parser.parse_args()

# Create a session with your profile
session = boto3.Session(profile_name='sudo')
# Now use this session to create your IAM client
iam = session.client('iam')

# Process each JSON file in the specified folder
for filename in os.listdir(args.folder):
    if filename.endswith('.json'):
        policy_name = os.path.splitext(filename)[0]
        file_path = os.path.join(args.folder, filename)        
        with open(file_path, 'r') as file:
            policy_document = json.load(file)
            create_policy(policy_name, policy_document)

            # Get the ARN of the newly created policy
            policy_arn = iam.get_policy(PolicyName=policy_name)['Policy']['Arn']

            # Attach the policy to all users
            users = iam.list_users()['Users']
            for user in users:
                attach_policy_to_user(policy_arn, user['UserName'])
    elif filename == 'allow-existing-policies.txt':
        # Process allow-existing-policies.txt file
        policy_list_file = os.path.join(args.folder, filename)
        if os.path.exists(policy_list_file):
            with open(policy_list_file, 'r') as file:
                for line in file:
                    policy_arn = line.strip()
                    if not policy_arn.startswith('arn:'):
                        policy_arn = f"arn:aws:iam::aws:policy/{policy_arn}"
                    # Attach the policy to all users
                    users = iam.list_users()['Users']
                    for user in users:
                        attach_policy_to_user(policy_arn, user['UserName'])        
