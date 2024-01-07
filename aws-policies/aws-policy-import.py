#! /usr/bin/env python3

import boto3, json, os, sys, argparse
from botocore.exceptions import ClientError

def create_policy(policy_name, policy_document):
    try:
        # Try to create the policy
        response = iam.create_policy(
            PolicyName=policy_name,
            PolicyDocument=json.dumps(policy_document)
        )
        policy_arn = response['Policy']['Arn']
        print(f"Policy {policy_name} created successfully.")
        return policy_arn
    
    except iam.exceptions.EntityAlreadyExistsException:

        # If the policy already exists, find its ARN and update it
        policies = iam.list_policies(Scope='Local')['Policies']
        policy_arn = next((p['Arn'] for p in policies if p['PolicyName'] == policy_name), None)

        policy_versions = iam.list_policy_versions(PolicyArn=policy_arn)['Versions']
        non_default_versions = [v for v in policy_versions if not v['IsDefaultVersion']]

        if non_default_versions:
            oldest_version = sorted(non_default_versions, key=lambda x: x['CreateDate'])[0]
            iam.delete_policy_version(PolicyArn=policy_arn, VersionId=oldest_version['VersionId'])
            print(f"Deleted oldest non-default policy version: {oldest_version['VersionId']}")

        if policy_arn:
            iam.create_policy_version(
                PolicyArn=policy_arn,
                PolicyDocument=json.dumps(policy_document),
                SetAsDefault=True
            )
            print(f"Policy {policy_name} already exists and was overwritten.")
            return policy_arn
        else:
            print(f"Failed to find ARN for existing policy {policy_name}")
            return None
    except ClientError as error:
        print(f"Error in creating/updating policy {policy_name}: {error}")
        sys.exit(1)


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
            policy_arn = create_policy(policy_name, policy_document)

            # Get the ARN of the newly created policy
            #policy_arn = iam.get_policy(PolicyName=policy_name)['Policy']['Arn']

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
