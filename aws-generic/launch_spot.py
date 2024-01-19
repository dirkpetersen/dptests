#! /usr/bin/env python3

import boto3
import time

class BigBadClass:
    def __init__(self):
        self.ec2 = boto3.client('ec2')

    def get_on_demand_price(self, instance_type, region):
        # Implement logic to fetch the on-demand price for the given instance type and region.
        # You might use AWS Price List API here.
        # For the sake of this example, I'm returning a hardcoded price.
        return 0.1  # Example price, replace with actual API call

    def launch_spot_instance(self, instance_type, ami_id, max_price, region):
        try:
            response = self.ec2.request_spot_instances(
                SpotPrice=str(max_price),
                InstanceCount=1,
                Type='one-time',
                LaunchSpecification={
                    'InstanceType': instance_type,
                    'ImageId': ami_id,
                }
            )
            return response
        except Exception as e:
            print(f"Failed to launch spot instance: {e}")
            return None

    def attempt_launch(self, instance_type, ami_id, region):
        on_demand_price = self.get_on_demand_price(instance_type, region)
        bid_price = on_demand_price / 3
        max_price = on_demand_price * 2 / 3

        for attempt in range(5):
            print(f"Trying to launch instance with bid price: {bid_price}")
            response = self.launch_spot_instance(instance_type, ami_id, bid_price, region)
            if response:
                print("Instance launched successfully.")
                return response
            else:
                bid_price += (max_price - bid_price) / 4
                if bid_price > max_price:
                    bid_price = max_price

            time.sleep(10)  # Wait before next attempt

        print("Unable to launch instance at desired price. Consider using on-demand.")
        return None

# Usage example
instance_type = "t3.micro"
ami_id = "ami-0b6e1e9474c3ee381"  # Replace with actual AMI ID
region = "us-west-2"
big_bad = BigBadClass()
big_bad.attempt_launch(instance_type, ami_id, region)

