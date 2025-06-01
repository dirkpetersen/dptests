import boto3
import json
import logging
from datetime import datetime

# Configure logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def lambda_handler(event, context):
    """
    Lambda function to delete unused EBS volumes.
    Runs nightly at 2:00 AM via EventBridge schedule.
    """
    
    ec2 = boto3.client('ec2')
    
    try:
        # Get all volumes
        response = ec2.describe_volumes()
        volumes = response['Volumes']
        
        deleted_volumes = []
        errors = []
        
        for volume in volumes:
            volume_id = volume['VolumeId']
            state = volume['State']
            
            # Check if volume is available (not attached)
            if state == 'available':
                try:
                    # Additional safety check - ensure no attachments
                    if not volume.get('Attachments', []):
                        logger.info(f"Deleting unused volume: {volume_id}")
                        
                        # Delete the volume
                        ec2.delete_volume(VolumeId=volume_id)
                        deleted_volumes.append(volume_id)
                        
                    else:
                        logger.info(f"Skipping volume {volume_id} - has attachments")
                        
                except Exception as e:
                    error_msg = f"Failed to delete volume {volume_id}: {str(e)}"
                    logger.error(error_msg)
                    errors.append(error_msg)
            else:
                logger.info(f"Skipping volume {volume_id} - state: {state}")
        
        # Log summary
        logger.info(f"Deleted {len(deleted_volumes)} unused volumes")
        if deleted_volumes:
            logger.info(f"Deleted volumes: {', '.join(deleted_volumes)}")
        
        if errors:
            logger.error(f"Errors encountered: {errors}")
        
        return {
            'statusCode': 200,
            'body': json.dumps({
                'message': f'Successfully processed {len(volumes)} volumes',
                'deleted_count': len(deleted_volumes),
                'deleted_volumes': deleted_volumes,
                'errors': errors,
                'timestamp': datetime.utcnow().isoformat()
            })
        }
        
    except Exception as e:
        error_msg = f"Lambda execution failed: {str(e)}"
        logger.error(error_msg)
        
        return {
            'statusCode': 500,
            'body': json.dumps({
                'error': error_msg,
                'timestamp': datetime.utcnow().isoformat()
            })
        }
