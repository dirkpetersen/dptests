#! /usr/bin/env python3

import boto3

core_account_id = '45488595XXXX'

def assume_role(account_id):
    sts_client = boto3.client('sts')
    assumed_role = sts_client.assume_role(
        RoleArn=f'arn:aws:iam::{account_id}:role/dns-creator-core',        
        RoleSessionName='CreateDNSRecordSession'
    )
    return assumed_role['Credentials']

def create_dns_record(account_id):
    credentials = assume_role(account_id)

    route53_client = boto3.client(
        'route53',
        aws_access_key_id=credentials['AccessKeyId'],
        aws_secret_access_key=credentials['SecretAccessKey'],
        aws_session_token=credentials['SessionToken'],
    )

    change_batch = {
        'Changes': [{
            'Action': 'CREATE',
            'ResourceRecordSet': {
                'Name': 'example3.r.internetxxxx.de',
                'Type': 'A',
                'TTL': 300,
                'ResourceRecords': [{'Value': '192.0.2.44'}]
            }
        }]
    }

    try:        
        hosted_zones = route53_client.list_hosted_zones()['HostedZones']
        if not hosted_zones:
            print('No hosted zones found, cannot create DNS record')
            return    
        # get the first found hosted zone id 
        hosted_zone_id = hosted_zones[0]['Id'].split('/')[-1] 
        print('Hosted Zone ID:', hosted_zone_id) 
        response = route53_client.change_resource_record_sets(
            HostedZoneId=hosted_zone_id,
            ChangeBatch=change_batch
        )
        print(response)        
    except Exception as e:    
        print(f'Error creating DNS record: {e}')
        
if __name__ == '__main__':   
    session = boto3.Session()
    print('My Account:', session.client('sts').get_caller_identity().get('Account'))    
    create_dns_record(core_account_id)

