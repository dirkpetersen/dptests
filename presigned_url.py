import datetime
import hashlib
import hmac
import urllib.parse

class AWSTranscribePresignedURL:
    def __init__(self, access_key, secret_key, session_token, region):
        self.access_key = access_key
        self.secret_key = secret_key
        self.session_token = session_token
        self.region = region
        self.service = 'transcribe'

    def get_request_url(self, sample_rate, language_code, media_encoding, 
                       number_of_channels=1, enable_channel_identification=False):
        """Generate a pre-signed URL for Amazon Transcribe WebSocket."""
        
        # Request parameters
        endpoint = f"wss://transcribestreaming.{self.region}.amazonaws.com:8443"
        host = f"transcribestreaming.{self.region}.amazonaws.com:8443"
        
        # Create date strings
        t = datetime.datetime.utcnow()
        amz_date = t.strftime('%Y%m%dT%H%M%SZ')
        date_stamp = t.strftime('%Y%m%d')
        
        # Task 1: Create a Canonical Request
        canonical_uri = "/stream-transcription-websocket"
        
        # Create the canonical query string
        credential_scope = f"{date_stamp}/{self.region}/{self.service}/aws4_request"
        
        canonical_querystring = {
            "X-Amz-Algorithm": "AWS4-HMAC-SHA256",
            "X-Amz-Credential": f"{self.access_key}/{credential_scope}",
            "X-Amz-Date": amz_date,
            "X-Amz-Expires": "300",
            "X-Amz-SignedHeaders": "host",
            "language-code": language_code,
            "media-encoding": media_encoding,
            "sample-rate": str(sample_rate)
        }
        
        if self.session_token:
            canonical_querystring["X-Amz-Security-Token"] = self.session_token
            
        if number_of_channels > 1:
            canonical_querystring["number-of-channels"] = str(number_of_channels)
            if enable_channel_identification:
                canonical_querystring["enable-channel-identification"] = "true"
        
        # Sort and encode query string
        canonical_querystring = "&".join([
            f"{urllib.parse.quote(k, safe='')}={urllib.parse.quote(str(v), safe='')}"
            for k, v in sorted(canonical_querystring.items())
        ])
        
        # Create the canonical headers
        canonical_headers = f"host:{host}\n"
        
        # Create the signed headers
        signed_headers = "host"
        
        # Create payload hash
        payload_hash = hashlib.sha256(b"").hexdigest()
        
        # Combine elements to create canonical request
        canonical_request = "\n".join([
            "GET",
            canonical_uri,
            canonical_querystring,
            canonical_headers,
            signed_headers,
            payload_hash
        ])
        
        # Task 2: Create the String to Sign
        algorithm = "AWS4-HMAC-SHA256"
        string_to_sign = "\n".join([
            algorithm,
            amz_date,
            credential_scope,
            hashlib.sha256(canonical_request.encode('utf-8')).hexdigest()
        ])
        
        # Task 3: Calculate the Signature
        def sign(key, msg):
            return hmac.new(key, msg.encode('utf-8'), hashlib.sha256).digest()
        
        k_date = sign(f"AWS4{self.secret_key}".encode('utf-8'), date_stamp)
        k_region = sign(k_date, self.region)
        k_service = sign(k_region, self.service)
        k_signing = sign(k_service, "aws4_request")
        signature = hmac.new(k_signing, string_to_sign.encode('utf-8'), hashlib.sha256).hexdigest()
        
        # Task 4: Add signature to query string
        canonical_querystring += f"&X-Amz-Signature={signature}"
        
        # Create request URL
        request_url = f"{endpoint}{canonical_uri}?{canonical_querystring}"
        
        return request_url
