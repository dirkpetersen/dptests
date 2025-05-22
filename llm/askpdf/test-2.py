#! /usr/bin/env python3

import boto3
import sys
from botocore.exceptions import ClientError

region = 'us-east-1'
print('Using region: ', region)

def document_conversation(bedrock_client,
                     model_id,
                     input_text,
                     input_document):
    """
    Sends a message to a model.
    Args:
        bedrock_client: The Boto3 Bedrock runtime client.
        model_id (str): The model ID to use.
        input text : The input message.
        input_document : The input document.

    Returns:
        response (JSON): The conversation that the model generated.

    """

    print(f"Generating message with model {model_id}")

    # Message to send.
    
    with open(input_document, "rb") as f:
        doc_bytes = f.read()

    message = {
        "role": "user",
        "content": [
            {
                "text": input_text
            },
            {
                "document": {
                    "name": "MyDocument",
                    "format": "pdf",
                    "source": {
                        "bytes": doc_bytes
                    }
                }
            }
        ]
    }

    messages = [message]

    # Send the message.
    response = bedrock_client.converse(
        modelId=model_id,
        messages=messages
    )

    return response


#model_id = "mistral.mistral-large-2402-v1:0"
model_id = "us.amazon.nova-pro-v1:0"

input_text = "What's in this document?"
input_document = 'ChrisM.pdf'

try:

    bedrock_client = boto3.client(service_name="bedrock-runtime", region_name=region)
    print(f"Bedrock client configured for region: {bedrock_client.meta.region_name}")

    response = document_conversation(
        bedrock_client, model_id, input_text, input_document)

    output_message = response['output']['message']

    print(f"Role: {output_message['role']}")

    for content in output_message['content']:
        print(f"Text: {content['text']}")

    token_usage = response['usage']
    print(f"Input tokens:  {token_usage['inputTokens']}")
    print(f"Output tokens:  {token_usage['outputTokens']}")
    print(f"Total tokens:  {token_usage['totalTokens']}")
    print(f"Stop reason: {response['stopReason']}")

except ClientError as err:
    message = err.response['Error']['Message']
    print(f"A client error occured: {message}")

else:
    print(
        f"Finished generating text with model {model_id}.")
