#! /usr/bin/env python3

import boto3, botocore
import json

# Initialize the Bedrock runtime client
client = boto3.client('bedrock-runtime', region_name='us-west-2')

# Your model's ARN
model_id = 'arn:aws:bedrock:us-west-2:405644541454:imported-model/d1jzq2drv2cb'
#model_id = 'arn:aws:bedrock:us-west-2:405644541454:imported-model/deepseek-r1-distill-llama-70b'

# Example inference call
def invoke_model(prompt):
    response = client.invoke_model(
        modelId=model_id,
        body=json.dumps({'prompt': prompt}),
        accept='application/json',
        contentType='application/json'
    )

    return json.loads(response['body'].read().decode('utf-8'))

# Example usage
try:
    result = invoke_model("Explain quantum computing in simple terms.")
    print(result)
except botocore.exceptions.ModelNotReadyException as e:
    print(f'Model {model_id} is not readly yet. Try again later. Error: {e}')
    
