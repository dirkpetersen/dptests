#! /usr/bin/env python3

import boto3, json

client_account_id = '70174744XXXX'
aws_profile = 'dipeXX'

def create_route53_policy(iam, hosted_zone, core_account):
    policy_document = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": [
                    "route53:ChangeResourceRecordSets"
                ],
                "Resource": f"arn:aws:route53:::hostedzone/{hosted_zone}",
                "Condition": {
                    "ForAllValues:StringEquals": {
                        "route53:ChangeResourceRecordSetsActions": [
                            "CREATE"
                        ]
                    }
                }
            },
            {
                "Effect": "Allow",
                "Action": [
                    "route53:ListResourceRecordSets"
                ],
                "Resource": f"arn:aws:route53:::hostedzone/{hosted_zone}"
            },
            {
                "Effect": "Allow",
                "Action": [
                    "route53:ListHostedZones"
                ],
                "Resource": "*"
            }                        
        ]
    }

    # policy_document = {
    #     "Version": "2012-10-17",
    #     "Statement": [
    #         {
    #             "Effect": "Allow",
    #             "Action": [
    #                 "route53:*"
    #             ],
    #              "Resource": "*"
    #         }
    #     ]
    # }    

    try:
        response = iam.create_policy(
            PolicyName='CrossAccountDNSChangePolicy',
            PolicyDocument=json.dumps(policy_document)
        )
    except Exception as e:
        print(f"Policy might already exists: {e}")
        response = iam.get_policy(PolicyArn=f"arn:aws:iam::{core_account}:policy/CrossAccountDNSChangePolicy")
        
    return response['Policy']['Arn']

def create_role_with_trust(iam, client_account, policy_arn):
    trust_policy = {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {
                    "AWS": [
                        f"arn:aws:iam::{client_account}:root" # allow all users in the account to use the role
                    ]
                },                
                "Action": "sts:AssumeRole"
            }
        ]
    }
    print(f'Creating role "dns-creator-core" for account {client_account} to assume the role')
    try:
        role = iam.create_role(
            RoleName='dns-creator-core',
            AssumeRolePolicyDocument=json.dumps(trust_policy)
        )
    except Exception as e:
        print(f"Role might already exists: {e}")
        role = iam.get_role(RoleName='dns-creator-core')

    print(f'Attaching policy CrossAccountDNSChangePolicy to role "dns-creator-core"')
    iam.attach_role_policy(
        RoleName='dns-creator-core',
        PolicyArn=policy_arn
    )

    return role['Role']['Arn']

def main():
    session = boto3.Session(profile_name=aws_profile)
    iam = session.client('iam')    
    r53 = session.client('route53')
    core_account = session.client('sts').get_caller_identity().get('Account')
    print('My Account:', core_account)

    hosted_zones = r53.list_hosted_zones()['HostedZones']
    if not hosted_zones:
        print('No hosted zones found, cannot create DNS record')
        return    
    hosted_zone = hosted_zones[0]['Id'].split('/')[-1]

    policy_arn = create_route53_policy(iam, hosted_zone, core_account)
    print(f"Policy ARN: {policy_arn}")
    role_arn = create_role_with_trust(iam, client_account_id, policy_arn)
    print(f"Role ARN: {role_arn}")
    

if __name__ == '__main__':
    main()

