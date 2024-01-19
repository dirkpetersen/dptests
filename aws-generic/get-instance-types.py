#! /usr/bin/env python3
import boto3

class BigBadClass:
    def list_7th_gen_instance_types(self, region='us-west-2'):
        ec2 = boto3.client('ec2', region_name=region)

        # Retrieve all instance types
        paginator = ec2.get_paginator('describe_instance_types')
        page_iterator = paginator.paginate()

        # Filter for 7th generation instance types
        gen_7_instance_types = []
        for page in page_iterator:
            for instance_type in page['InstanceTypes']:
                if instance_type['InstanceType'].startswith(('m7', 'c7', 'r7')):
                    gen_7_instance_types.append(instance_type['InstanceType'])

        return gen_7_instance_types

# Example usage
big_bad = BigBadClass()
gen_7_instance_types = big_bad.list_7th_gen_instance_types()
print(gen_7_instance_types)


