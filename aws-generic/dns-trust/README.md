
# DNS Route53 cross account setup 

## dns-core-permissions.py

- sets up a role and a policy with access to a DNS Hosted Zone in a central 'core' account
- you need to set 'client_account_id' and 'aws_profile' variables in the code 

## dns-client-permissions.py

- allows members of a group to access to dns role in the 'core' account from the client account 
- you need to set 'core_account_id', 'aws_profile' and 'group_name' variables

## dns-record-as-client.py

- creates a new DNS A record in first found zone id. CREATE is allowed, UPSERT NOT
- you need to set 'core_account_id' variable

