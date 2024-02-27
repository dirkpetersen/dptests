#! /usr/bin/env python3

import boto3, json

core_account_id = "45488595XXXX"
aws_profile = 'sudo'
group_name = "dpcriDevs"

# THIS IS ONLY FOR ROLE CHAINING
# 
# def create_cross_account_role(iam, server_account):
#     # Initialize a boto3 client for IAM
    
#     # Define the trust policy for the dns-creator role to be assumed by peterdir
#     trust_policy = {
#         "Version": "2012-10-17",
#         "Statement": [
#             {
#                 "Effect": "Allow",
#                 "Principal": {"AWS": "arn:aws:iam::701747442027:user/peterdir"},
#                 "Action": "sts:AssumeRole"
#             }
#         ]
#     }

#     # Define the policy allowing dns-creator to assume dns-creator-core in another account
#     assume_role_policy = {
#         "Version": "2012-10-17",
#         "Statement": [
#             {
#                 "Effect": "Allow",
#                 "Action": "sts:AssumeRole",
#                 "Resource": f"arn:aws:iam::{server_account}:role/dns-creator-core"
#             }
#         ]
#     }

#     # Create the dns-creator role with the trust policy
#     role = iam.create_role(
#         RoleName='dns-creator',
#         AssumeRolePolicyDocument=json.dumps(trust_policy)
#     )

#     # Attach the policy allowing dns-creator to assume dns-creator-core
#     iam.put_role_policy(
#         RoleName='dns-creator',
#         PolicyName='AssumeDnsCreatorSrvPolicy',
#         PolicyDocument=json.dumps(assume_role_policy)
#     )

#     return role['Policy']['Arn']


def create_remote_role_policy(iam, core_account_id):
    policy_document = {        
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": "sts:AssumeRole",
                "Resource": f"arn:aws:iam::{core_account_id}:role/dns-creator-core"
            }
        ]
    }
    try:
        response = iam.create_policy(
            PolicyName='RemoteRolePolicy',
            PolicyDocument=json.dumps(policy_document)
        )    
        newarn = response['Policy']['Arn']
        print(f'Created policy {newarn}')

    except Exception as e:
        print(f"Policy might already exists: {e}")       
        newarn = None

    return newarn

# def attach_policy_to_user(iam, user, policy_arn):
#     result = iam.attach_user_policy(
#         UserName=user,
#         PolicyArn=policy_arn
#     )
#     print(f'Attached policy {policy_arn} to user {user}')
#     return result

def attach_policy_to_group(iam, group, policy_arn):

    try:
        result = iam.attach_group_policy(
            GroupName=group,
            PolicyArn=policy_arn
        )
        print(f'Attached policy {policy_arn} to group {group}')
    except Exception as e:      
        print(f'Error attaching policy {policy_arn} to group {group}: {e}')
        result = None
    return result

def main():

    session = boto3.Session(profile_name=aws_profile)
    iam = session.client('iam')    
    print('My Account:', session.client('sts').get_caller_identity().get('Account'))

    #role_arn = create_cross_account_role(iam, account_id)
    pol_arn = create_remote_role_policy(iam, core_account_id)    
    print(f"Policy ARN: {pol_arn}")
    if pol_arn:
        #ret = attach_policy_to_user(iam, 'peterdir', pol_arn)
        ret = attach_policy_to_group(iam, group_name, pol_arn)
    else:
        print('no new policy, cannot run attach_policy_to_group')
    
if __name__ == '__main__':
    main()
