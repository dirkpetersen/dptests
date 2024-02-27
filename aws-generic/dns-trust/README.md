
# DNS Route53 cross account setup 

## dns-core-permissions.py

sets up a role and a policy with access to a DNS Hosted Zone in a central 'core' account 

## dns-client-permissions.py

allows members of a group to access to dns role in the 'core' account from the client account 

## dns-record-as-client.py

creates a new DNS A record in first found zone id. CREATE is allowed, UPSERT NOT 

