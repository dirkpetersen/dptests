#! /usr/bin/env python3

import boto3, datetime
from operator import itemgetter

class BigBadClass:

    def __init__(self):
        self.cpu_types = {
            "graviton-2": ('m6g','c6g', 'c6gn', 't4g' ,'g5g'),
            "graviton-3": ('m7g', 'c7g', 'c7gn'),
            "epyc-gen-1": ('t3a'),
            "epyc-gen-2": ('c5a', 'm5a', 'r5a', 'g4ad', 'p4', 'inf2', 'g5'),
            "epyc-gen-3": ('m6a', 'c6a', 'r6a', 'p5'),
            "epyc-gen-4": ('c7a', 'm7a', 'r7a'),
            "xeon-gen-1": ('m4', 'c4', 't2', 'r4', 'p3' ,'p2', 'f1', 'g3', 'i3en'),
            "xeon-gen-2": ('m5', 'c5', 'c5n', 'm5n', 'm5zn', 'r5', 't3', 't3n', 'dl1', 'inf1', 'g4dn', 'vt1'),
            "xeon-gen-3": ('m6i', 'c6i', 'm6in', 'c6in', 'r6i', 'r6id', 'r6idn', 'r6in', 'trn1'),
            "xeon-gen-4": ('c7i', 'm7i', 'm7i-flex'),
            "core-i7-mac": ('mac1')
        }

    def _ec2_describe_instance_families(self, cpu_type, vcpus=1, memory_gb=1, region='us-west-2'):

        ec2 = boto3.client('ec2', region_name=region)
        instance_families = self.cpu_types[cpu_type]

        filtered_instance_families = []

        # Retrieve all instance types
        paginator = ec2.get_paginator('describe_instance_types')
        page_iterator = paginator.paginate()

        for page in page_iterator:
            for instance_type in page['InstanceTypes']:
                if instance_type['InstanceType'].startswith(instance_families) and \
                   instance_type['VCpuInfo']['DefaultVCpus'] >= vcpus and \
                   instance_type['MemoryInfo']['SizeInMiB'] >= memory_gb * 1024:
                    filtered_instance_families.append(instance_type)

        #instance_ids = [i['InstanceType'] for i in filtered_instance_families]
        #print('\nfiltered_instance_families:', instance_ids)

        return filtered_instance_families

    def _ec2_get_cheapest_spot_instance(self, cpu_type, vcpus=1, memory_gb=1, regions=['us-west-2', 'us-west-1', 'us-east-2', 'us-east-1', 'ca-central-1', 'eu-north-1', 'eu-central-1']):
        
        # Validate CPU type
        if cpu_type not in self.cpu_types:
            return "Invalid CPU type.", None, None    
      
        # Filter instances by vCPUs, memory, and CPU type
        filtered_instances = self._ec2_describe_instance_families(cpu_type, vcpus, memory_gb)
        if not filtered_instances:
            return "No instances match the criteria."
        
        instance_ids = [i['InstanceType'] for i in filtered_instances]
        print('\ninstance_ids:', instance_ids)
        start_time = datetime.datetime.utcnow() - datetime.timedelta(minutes=15)         
        cheapest_instance = None
        lowest_price = float('inf')                       

        for region in regions:
            print(f'Checking AWS region: {region}')
            ec2 = boto3.client('ec2', region_name=region)
            # Get current spot prices for filtered instances            
            spot_prices = ec2.describe_spot_price_history(
                StartTime=start_time,
                InstanceTypes=instance_ids,
                ProductDescriptions=['Linux/UNIX'],
                MaxResults=len(instance_ids)
            )
            for price in spot_prices['SpotPriceHistory']:
                if float(price['SpotPrice']) < lowest_price:
                    lowest_price = float(price['SpotPrice'])
                    cheapest_instance = price
            if cheapest_instance:
                print(f"   Cheapest Instance: {cheapest_instance['InstanceType']}, Price: {cheapest_instance['SpotPrice']}, Region: {cheapest_instance['AvailabilityZone']}")


        if cheapest_instance:
            print(f"*** Cheapest Instance: {cheapest_instance['InstanceType']}, Price: {cheapest_instance['SpotPrice']}, Region: {cheapest_instance['AvailabilityZone']}")
        else:
            print("No cheapest instance found.")            

        #print('\nspot_prices:', spot_prices)

        # Find the cheapest instance
        #cheapest_instance = min(spot_prices['SpotPriceHistory'], key=itemgetter('SpotPrice'))
        return cheapest_instance['InstanceType'], cheapest_instance['AvailabilityZone'], cheapest_instance['SpotPrice']

# Example usage
big_bad = BigBadClass()
vcpus = 4 # Number of vCPUs
memory_gb = 8  # Memory in GB
#cpu_type = 'xeon-gen-1'  # CPU type
#cpu_type = 'xeon-gen-2'
#cpu_type = 'xeon-gen-3'
#cpu_type = 'xeon-gen-4'
#cpu_type = 'graviton-2'  # CPU type
cpu_type = 'graviton-3'  # CPU type
#cpu_type = 'epyc-gen-1'
#cpu_type = 'epyc-gen-2'
#cpu_type = 'epyc-gen-3'
#cpu_type = 'epyc-gen-4'


#filtered_instances = big_bad.describe_instance_families(cpu_type, 32)
# Display some information about the filtered instances
#for instance in filtered_instances:
#    print(instance['InstanceType'])

cheapest_instance = big_bad._ec2_get_cheapest_spot_instance(cpu_type, vcpus, memory_gb)
print(f"Cheapest instance: {cheapest_instance}")

