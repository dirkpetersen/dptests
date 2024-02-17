
# Let's encrypt with AWS

## DNS first 

### create recordset.json

```
{
  "Comment": "Creating an A record",
  "Changes": [
    {
      "Action": "CREATE",
      "ResourceRecordSet": {
        "Name": "aws-eb.aws.internetchen.de",
        "Type": "A",
        "TTL": 300,
        "ResourceRecords": [
          {
            "Value": "35.xxx.xxx.xxx"
          }
        ]
      }
    }
  ]
}
```

register dns to first found hosted zone id with Route 53 and recordset.json 

```
first_hosted_zone_id=$(aws route53 list-hosted-zones | jq -r '.HostedZones[0].Id' | cut -d'/' -f3)
if [[ -n ${first_hosted_zone_id} ]]; then 
  aws route53 change-resource-record-sets --hosted-zone-id ${first_hosted_zone_id} --change-batch file://recordset.json
fi

```

you do not require a json file but can pass the json structure as a string

```
JSON='{"Changes":[{"Action":"CREATE","ResourceRecordSet":{"Name":"example.com","Type":"A","TTL":300,"ResourceRecords":[{"Value":"192.0.2.1"}]}}]}'
aws route53 change-resource-record-sets --hosted-zone-id ZONEID --change-batch $JSON
```

## Certbot

### register cert 
```
if ! sudo test -d /root/.aws; then
  # temp access to AWS creds for root user if not set for root
  sudo ln -s /home/ec2-user/.aws /root/.aws
fi
python3 -m venv le
. le/bin/activate
pip install certbot-dns-route53
sudo /home/ec2-user/le/bin/certbot certonly --dns-route53 --register-unsafely-without-email --agree-tos -d aws-eb.aws.internetchen.de
sudo rm -f /root/.aws
```

### renew cert 

```
if ! sudo test -d /root/.aws; then
  # temp access to AWS creds for root user if not set for root
  sudo ln -s /home/ec2-user/.aws /root/.aws
fi
sudo /home/ec2-user/le/bin/certbot renew --force-renewal --dns-route53
rm -f /root/.aws
```

### create systemd timer to renew every 75 days

```
# Create a systemd service file
cat <<EOT >> /etc/systemd/system/certbot-renew.service
[Unit]
Description=Let's Encrypt renewal

[Service]
Type=oneshot
ExecStart=/home/ec2-user/le/bin/certbot renew --force-renewal --dns-route53
EOT

# Create a systemd timer file
cat <<EOT >> /etc/systemd/system/certbot-renew.timer
[Unit]
Description=Run Let's Encrypt renewal every 75 days

[Timer]
OnBootSec=1h
OnUnitActiveSec=75d

[Install]
WantedBy=timers.target
EOT

# Enable and start the timer
systemctl enable certbot-renew.timer
systemctl start certbot-renew.timer
```

### Problems 

```
Requesting a certificate for aws-eb.aws.internetchen.de
An unexpected error occurred:
Error creating new order :: too many certificates (5) already issued for this exact set of domains in the last 168 hours: aws-eb.aws.internetchen.de, retry after 2024-02-18T02:58:58Z: see https://letsencrypt.org/docs/duplicate-certificate-limit/
Ask for help or search for solutions at https://community.letsencrypt.org. See the logfile /var/log/letsencrypt/letsencrypt.log or re-run Certbot with -v for more details.
```


