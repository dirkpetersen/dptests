# DNS Route53 cross account setup 

dns_creator_core=
core_account_id=12345678912
client_account_id=

## shared parameters

we do not want to document / store the account id of the core account anywhere, the core account
can share account id or the role arn of the dns-creator role: 
https://docs.aws.amazon.com/systems-manager/latest/userguide/parameter-store-shared-parameters.html

## dns-core-permissions.py

- sets up a role and a policy with access to a DNS Hosted Zone in a central 'core' account
- you need to set 'client_account_id' and 'aws_profile' variables in the code 

## dns-client-permissions.py

- allows members of a group to access to dns role in the 'core' account from the client account 
- you need to set 'core_account_id', 'aws_profile' and 'group_name' variables

## dns-record.py

- creates a new DNS A record in first found zone id. CREATE is allowed, UPSERT NOT
- you need to set 'core_account_id' variable

