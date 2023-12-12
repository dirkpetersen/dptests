#! /usr/bin/env python3

import boto3
import pandas as pd
from datetime import datetime
from datetime import timedelta
ec2c = boto3.client('ec2')
ec2r = boto3.resource('ec2')

#### The rest of this code maps the instance details to spot price in case you are looking for certain memory or cpu
paginator = ec2c.get_paginator('describe_instance_types')
response_iterator = paginator.paginate( )

df_hold_list = []
for page in response_iterator:
    df_hold_list.append(pd.DataFrame(page['InstanceTypes']))

df_instance_specs = pd.concat(df_hold_list, axis=0).reset_index(drop=True)
df_instance_specs['Spot'] = df_instance_specs['SupportedUsageClasses'].apply(lambda x: 1 if 'spot' in x else 0)
df_instance_spot_specs = df_instance_specs.loc[df_instance_specs['Spot']==1].reset_index(drop=True)

#unapck memory and cpu dictionaries
df_instance_spot_specs['MemSize'] = df_instance_spot_specs['MemoryInfo'].apply(lambda x: x.get('SizeInMiB'))
df_instance_spot_specs['vCPUs'] = df_instance_spot_specs['VCpuInfo'].apply(lambda x: x.get('DefaultVCpus'))
df_instance_spot_specs['Processor'] = df_instance_spot_specs['ProcessorInfo'].apply(lambda x: x.get('SupportedArchitectures'))

#look at instances only between 30MB and 70MB
instance_list = df_instance_spot_specs['InstanceType'].unique().tolist()

#---------------------------------------------------------------------------------------------------------------------
# You can use this section by itself to get the instancce type and availability zone and loop through the instance you want
# just modify instance_list with one instance you want informatin for
#look only in us-east-1
client = boto3.client('ec2', region_name='eu-central-1')
prices = client.describe_spot_price_history(
    InstanceTypes=instance_list,
    ProductDescriptions=['Linux/UNIX', 'Linux/UNIX (Amazon VPC)'],
    StartTime=(datetime.now() -
               timedelta(hours=1)).isoformat(),
               # AvailabilityZone='us-east-1a'
    MaxResults=1000)

df_spot_prices = pd.DataFrame(prices['SpotPriceHistory'])
df_spot_prices['SpotPrice'] = df_spot_prices['SpotPrice'].astype('float')
df_spot_prices.sort_values('SpotPrice', inplace=True)
#---------------------------------------------------------------------------------------------------------------------

# merge memory size and cpu information into this dataframe
df_spot_instance_options = df_spot_prices[['AvailabilityZone', 'InstanceType', 'SpotPrice']].merge(df_instance_spot_specs[['InstanceType', 'MemSize', 'vCPUs',
                                            'CurrentGeneration', 'Processor']], left_on='InstanceType', right_on='InstanceType')


print(df_spot_instance_options)
