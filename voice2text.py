import sys
import os
import time
import threading
import uuid
import boto3
import pystray
import requests
from PIL import Image
import pygetwindow as gw
import pyautogui
import sounddevice as sd
import soundfile as sf
import tempfile
from mic import record_audio

class TranscriptionApp:
    def __init__(self):
        self.recording = False
        self.client = boto3.client('transcribe')
        self.s3_client = boto3.client('s3')
        self.bucket_name = 'voice-transcribe-temp'  # You need to create this bucket
        self.setup_tray()
        self.ensure_bucket_exists()
    
    def ensure_bucket_exists(self):
        try:
            self.s3_client.head_bucket(Bucket=self.bucket_name)
        except:
            self.s3_client.create_bucket(
                Bucket=self.bucket_name,
                CreateBucketConfiguration={'LocationConstraint': self.s3_client.meta.region_name}
            )
        
    def setup_tray(self):
        # Create a simple icon (you may want to replace this with a proper icon file)
        icon_image = Image.new('RGB', (64, 64), color='red')
        self.icon = pystray.Icon(
            "Transcriber",
            icon_image,
            "Voice Transcriber",
            menu=self.create_menu()
        )

    def create_menu(self):
        return pystray.Menu(
            pystray.MenuItem(
                "Start Recording",
                self.toggle_recording,
                checked=lambda item: self.recording
            ),
            pystray.MenuItem("Quit", self.quit_app)
        )

    def toggle_recording(self):
        self.recording = not self.recording
        if self.recording:
            threading.Thread(target=self.record_and_transcribe, daemon=True).start()
        
    def quit_app(self):
        self.recording = False
        self.icon.stop()

    def record_and_transcribe(self):
        while self.recording:
            # Create a temporary file for the audio
            with tempfile.NamedTemporaryFile(suffix='.wav', delete=False) as temp_file:
                temp_filename = temp_file.name
            
            # Record a 5-second chunk
            record_audio(5, filename=temp_filename)
            
            # Upload to S3
            s3_key = f"audio_{uuid.uuid4()}.wav"
            self.s3_client.upload_file(temp_filename, self.bucket_name, s3_key)
            s3_uri = f"s3://{self.bucket_name}/{s3_key}"

            # Start transcription job
            job_name = f"transcription_{uuid.uuid4()}"
            self.client.start_transcription_job(
                TranscriptionJobName=job_name,
                Media={'MediaFileUri': s3_uri},
                MediaFormat='wav',
                LanguageCode='en-US'
            )
            
            # Wait for transcription to complete
            while True:
                status = self.client.get_transcription_job(TranscriptionJobName=job_name)
                if status['TranscriptionJob']['TranscriptionJobStatus'] in ['COMPLETED', 'FAILED']:
                    break
                time.sleep(0.5)
            
            try:
                # Get transcription results
                if status['TranscriptionJob']['TranscriptionJobStatus'] == 'COMPLETED':
                    transcript_uri = status['TranscriptionJob']['Transcript']['TranscriptFileUri']
                    import requests
                    transcript_response = requests.get(transcript_uri)
                    transcript_data = transcript_response.json()
                    text = transcript_data['results']['transcripts'][0]['transcript']
                    
                    # Type the transcribed text into the active window
                    active_window = gw.getActiveWindow()
                    if active_window:
                        pyautogui.write(text)
            finally:
                # Cleanup
                os.unlink(temp_filename)
                self.s3_client.delete_object(Bucket=self.bucket_name, Key=s3_key)
            
            # Small delay before next recording
            time.sleep(0.1)

def main():
    app = TranscriptionApp()
    app.icon.run()

if __name__ == '__main__':
    main()
