#!/bin/bash

# Deploy Lambda function to delete unused EBS volumes
# Usage: ./deploy-lambda.sh

set -e

FUNCTION_NAME="delete-unused-ebs-volumes"
ROLE_NAME="DeleteUnusedVolumesRole"
POLICY_NAME="DeleteUnusedVolumesPolicy"
RULE_NAME="delete-unused-volumes-schedule"

echo "Starting deployment of Lambda function: $FUNCTION_NAME"

# Create IAM role if it doesn't exist
echo "Creating IAM role..."
aws iam create-role \
    --role-name $ROLE_NAME \
    --assume-role-policy-document '{
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Principal": {
                    "Service": "lambda.amazonaws.com"
                },
                "Action": "sts:AssumeRole"
            }
        ]
    }' 2>/dev/null || echo "Role already exists"

# Attach basic Lambda execution role
aws iam attach-role-policy \
    --role-name $ROLE_NAME \
    --policy-arn arn:aws:iam::aws:policy/service-role/AWSLambdaBasicExecutionRole

# Create custom policy for EC2 volume operations
echo "Creating custom IAM policy..."
aws iam create-policy \
    --policy-name $POLICY_NAME \
    --policy-document '{
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": [
                    "ec2:DescribeVolumes",
                    "ec2:DeleteVolume"
                ],
                "Resource": "*"
            }
        ]
    }' 2>/dev/null || echo "Policy already exists"

# Attach custom policy to role
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
aws iam attach-role-policy \
    --role-name $ROLE_NAME \
    --policy-arn arn:aws:iam::$ACCOUNT_ID:policy/$POLICY_NAME

# Wait for role to be ready
echo "Waiting for IAM role to be ready..."
sleep 10

# Create deployment package
echo "Creating deployment package..."
cd lambda
zip -r ../lambda-deployment.zip . -x "*.pyc" "__pycache__/*"
cd ..

# Create or update Lambda function
echo "Deploying Lambda function..."
aws lambda create-function \
    --function-name $FUNCTION_NAME \
    --runtime python3.9 \
    --role arn:aws:iam::$ACCOUNT_ID:role/$ROLE_NAME \
    --handler delete_unused_volumes.lambda_handler \
    --zip-file fileb://lambda-deployment.zip \
    --description "Deletes unused EBS volumes nightly at 2:00 AM" \
    --timeout 300 \
    --memory-size 128 2>/dev/null || \
aws lambda update-function-code \
    --function-name $FUNCTION_NAME \
    --zip-file fileb://lambda-deployment.zip

# Create EventBridge rule for nightly execution at 2:00 AM UTC
echo "Creating EventBridge schedule rule..."
aws events put-rule \
    --name $RULE_NAME \
    --schedule-expression "cron(0 2 * * ? *)" \
    --description "Trigger delete unused volumes Lambda at 2:00 AM daily"

# Add permission for EventBridge to invoke Lambda
echo "Adding EventBridge permission to Lambda..."
aws lambda add-permission \
    --function-name $FUNCTION_NAME \
    --statement-id "AllowExecutionFromEventBridge" \
    --action "lambda:InvokeFunction" \
    --principal events.amazonaws.com \
    --source-arn arn:aws:events:$(aws configure get region):$ACCOUNT_ID:rule/$RULE_NAME 2>/dev/null || echo "Permission already exists"

# Add Lambda as target to EventBridge rule
echo "Adding Lambda as target to EventBridge rule..."
aws events put-targets \
    --rule $RULE_NAME \
    --targets "Id"="1","Arn"="arn:aws:lambda:$(aws configure get region):$ACCOUNT_ID:function:$FUNCTION_NAME"

# Clean up deployment package
rm -f lambda-deployment.zip

echo "Deployment completed successfully!"
echo "Function Name: $FUNCTION_NAME"
echo "Schedule: Daily at 2:00 AM UTC"
echo "Next run: $(date -d 'tomorrow 02:00' '+%Y-%m-%d %H:%M:%S UTC')"
echo ""
echo "To test the function manually, run:"
echo "aws lambda invoke --function-name $FUNCTION_NAME output.json && cat output.json"
