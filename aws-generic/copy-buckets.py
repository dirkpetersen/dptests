
import boto3
import concurrent.futures
from botocore.exceptions import ClientError

def copy_object(s3, source_bucket, destination_bucket, obj):
    try:
        # Check if the object exists in the destination bucket
        dest_obj = s3.head_object(Bucket=destination_bucket, Key=obj['Key'], RequestPayer='requester')
        # Compare ETags (remove quotation marks from ETags if necessary)
        if dest_obj['ETag'] == obj['ETag']:
            print(f"Skipping {obj['Key']}, already exists in destination bucket")
            return
    except ClientError:
        # Object does not exist in the destination bucket
        pass

    # Copy object with Requester Pays option
    copy_source = {'Bucket': source_bucket, 'Key': obj['Key']}
    s3.copy(copy_source, destination_bucket, obj['Key'], 
            ExtraArgs={'RequestPayer': 'requester', 'StorageClass': tier})
    print(f"Copied {obj['Key']} from {source_bucket} to {destination_bucket}")

def copy_s3_objects(source_bucket, destination_bucket, max_workers=100):
    s3 = boto3.client('s3')

    paginator = s3.get_paginator('list_objects_v2')
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Iterate over each page in the paginator
        for page in paginator.paginate(Bucket=source_bucket, RequestPayer='requester'):
            if 'Contents' in page:
                # Submit each object copy to the thread pool
                futures = [executor.submit(copy_object, s3, source_bucket, destination_bucket, obj) 
                           for obj in page['Contents']]
                # Wait for all submitted futures to complete
                concurrent.futures.wait(futures)

# Replace with your source and destination bucket names
source_bucket_name = 'easybuild-cache-ohsu'
destination_bucket_name = 'easybuild-cache'

copy_s3_objects(source_bucket_name, destination_bucket_name)

