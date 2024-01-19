#! /usr/bin/env python3

import boto3
import json

class BigBadClass:
    def __init__(self):
        self.pricing_client = boto3.client('pricing', region_name='us-east-1')

    def get_on_demand_price(self, instance_type, region):
        region_map = {
            'us-east-1': 'US East (N. Virginia)',
            'us-west-1': 'US West (N. California)',
            'us-west-2': 'US West (Oregon)',
            # ... Add other regions as needed
        }
        try:
            response = self.pricing_client.get_products(
                ServiceCode='AmazonEC2',
                Filters=[
                    {'Type': 'TERM_MATCH', 'Field': 'instanceType', 'Value': instance_type},
                    {'Type': 'TERM_MATCH', 'Field': 'location', 'Value': region_map.get(region, '')},
                    {'Type': 'TERM_MATCH', 'Field': 'preInstalledSw', 'Value': 'NA'},
                    {'Type': 'TERM_MATCH', 'Field': 'operatingSystem', 'Value': 'Linux'},
                    {'Type': 'TERM_MATCH', 'Field': 'tenancy', 'Value': 'Shared'},
                    {'Type': 'TERM_MATCH', 'Field': 'capacitystatus', 'Value': 'Used'},
                ],
                MaxResults=1
            )
            price_list = [json.loads(price_str) for price_str in response['PriceList']]
            on_demand_price = float(price_list[0]['terms']['OnDemand'][list(price_list[0]['terms']['OnDemand'])[0]]['priceDimensions'][list(price_list[0]['terms']['OnDemand'][list(price_list[0]['terms']['OnDemand'])[0]]['priceDimensions'])[0]]['pricePerUnit']['USD'])
            return on_demand_price
        except Exception as e:
            print(f"Error getting on-demand price: {e}")
            return None

    @staticmethod
    def region_map(aws_region):
        region_map = {
            'us-east-1': 'US East (N. Virginia)',
            'us-west-1': 'US West (N. California)',
            'us-west-2': 'US West (Oregon)',
            # ... Add other regions as needed
        }
        return region_map.get(aws_region, '')

    # ... Other methods from the previous example

# Example usage
instance_type = "t4g.nano"
region = "us-west-2"
big_bad = BigBadClass()
price = big_bad.get_on_demand_price(instance_type, region)
print(f"On-demand price for {instance_type} in {region}: ${price}")

