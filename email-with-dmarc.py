#! /usr/bin/env python3

import boto3
from botocore.exceptions import ClientError
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import json

class SESEmailSender:
    def __init__(self, region_name='us-west-2'):
        """
        Initialize the SES client
        
        Args:
            region_name (str): AWS region name where SES is configured
        """
        self.client = boto3.client('ses', region_name=region_name)
        
    def create_message(self, sender, recipient, subject, body_text, body_html=None):
        """
        Create a MIME message that complies with DMARC
        
        Args:
            sender (str): Email address of the sender
            recipient (str): Email address of the recipient
            subject (str): Subject of the email
            body_text (str): Plain text version of the email body
            body_html (str): HTML version of the email body (optional)
        
        Returns:
            MIMEMultipart: The constructed MIME message
        """
        msg = MIMEMultipart('alternative')
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = recipient
        
        # Always include a plain text version
        part1 = MIMEText(body_text, 'plain')
        msg.attach(part1)
        
        # Optionally include an HTML version
        if body_html:
            part2 = MIMEText(body_html, 'html')
            msg.attach(part2)
            
        return msg
        
    def send_email(self, sender, recipient, subject, body_text, body_html=None):
        """
        Send an email using Amazon SES
        
        Args:
            sender (str): Email address of the sender
            recipient (str): Email address of the recipient
            subject (str): Subject of the email
            body_text (str): Plain text version of the email body
            body_html (str): HTML version of the email body (optional)
            
        Returns:
            dict: The response from the SES API
        
        Raises:
            ClientError: If there's an error sending the email
        """
        try:
            # Create the MIME message
            msg = self.create_message(sender, recipient, subject, body_text, body_html)
            
            # Convert the message to a raw format that SES can send
            raw_message = {
                'Data': msg.as_string()
            }
            
            # Send the email
            response = self.client.send_raw_email(
                Source=sender,
                Destinations=[recipient],
                RawMessage=raw_message
            )
            
            return response
            
        except ClientError as e:
            print(f"Error sending email: {e.response['Error']['Message']}")
            raise

# Example usage
if __name__ == "__main__":
    # Initialize the sender
    email_sender = SESEmailSender()
    
    # Example email content
    sender_email = "peterdir+ses@oregonstate.edu"
    recipient_email = "dipeit@gmail.com"
    subject = "Test Email with DMARC Compliance"
    
    text_content = """
    This is a test email sent using Amazon SES.
    This is the plain text version.
    """
    
    html_content = """
    <html>
    <head></head>
    <body>
        <h1>Test Email</h1>
        <p>This is a test email sent using Amazon SES.</p>
        <p>This is the HTML version.</p>
    </body>
    </html>
    """
    
    try:
        # Send the email
        response = email_sender.send_email(
            sender_email,
            recipient_email,
            subject,
            text_content,
            html_content
        )
        print(f"Email sent! Message ID: {response['MessageId']}")
        
    except ClientError as e:
        print(f"Failed to send email: {e}")
