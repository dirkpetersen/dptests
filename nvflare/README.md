# NVflare federated learning on AWS


Links: 

https://catalog.workshops.aws/nvflareonaws/en-US/deployment





Launch Single Account: 

catalog.workshops.aws/nvflareonaws/en-US/deployment


Single Account 
https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=nvflarestack&templateURL=https://monailabelonaws.s3.amazonaws.com/flareonaccountall.yaml


Multi-Account Network

https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=vpcstack&templateURL=https://monailabelonaws.s3.amazonaws.com/master_vpc.yml


Network Client 1:
https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=vpcstack&templateURL=https://monailabelonaws.s3.amazonaws.com/client1_vpc.yml

Network Client 2:
https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=vpcstack&templateURL=https://monailabelonaws.s3.amazonaws.com/client2_vpc.yml

Master Instance
https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=instancestack&templateURL=https://monailabelonaws.s3.amazonaws.com/master_instance.yml


Client Instance:
https://console.aws.amazon.com/cloudformation/home?region=us-east-1#/stacks/new?stackName=instancestack&templateURL=https://monailabelonaws.s3.amazonaws.com/client_instance.yml


