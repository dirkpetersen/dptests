import sys
import os
import asyncio
import uuid
import boto3
import pystray
import requests
import queue
import aiohttp
import json
import threading
import time
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
        self.audio_queue = queue.Queue()
        self.setup_tray()
        self.ensure_bucket_exists()
        self.loop = asyncio.new_event_loop()
    
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
            threading.Thread(target=self.start_async_loop, daemon=True).start()
            threading.Thread(target=self.record_audio_chunks, daemon=True).start()
        
    def quit_app(self):
        self.recording = False
        self.icon.stop()

    def start_async_loop(self):
        asyncio.set_event_loop(self.loop)
        self.loop.run_until_complete(self.process_audio_queue())

    def record_audio_chunks(self):
        while self.recording:
            with tempfile.NamedTemporaryFile(suffix='.wav', delete=False) as temp_file:
                temp_filename = temp_file.name
            
            # Record a 3-second chunk
            record_audio(3, filename=temp_filename)
            self.audio_queue.put(temp_filename)
            
            # Small overlap to ensure continuous recording
            time.sleep(2.8)

    async def process_audio_chunk(self, temp_filename):
        try:
            # Upload to S3
            s3_key = f"audio_{uuid.uuid4()}.wav"
            await self.loop.run_in_executor(
                None, 
                self.s3_client.upload_file,
                temp_filename, 
                self.bucket_name, 
                s3_key
            )
            s3_uri = f"s3://{self.bucket_name}/{s3_key}"

            # Start transcription job
            job_name = f"transcription_{uuid.uuid4()}"
            await self.loop.run_in_executor(
                None,
                lambda: self.client.start_transcription_job(
                    TranscriptionJobName=job_name,
                    Media={'MediaFileUri': s3_uri},
                    MediaFormat='wav',
                    LanguageCode='en-US'
                )
            )
            
            # Wait for transcription to complete
            while True:
                status = await self.loop.run_in_executor(
                    None,
                    lambda: self.client.get_transcription_job(TranscriptionJobName=job_name)
                )
                if status['TranscriptionJob']['TranscriptionJobStatus'] in ['COMPLETED', 'FAILED']:
                    break
                await asyncio.sleep(0.5)
            
            if status['TranscriptionJob']['TranscriptionJobStatus'] == 'COMPLETED':
                transcript_uri = status['TranscriptionJob']['Transcript']['TranscriptFileUri']
                async with aiohttp.ClientSession() as session:
                    async with session.get(transcript_uri) as response:
                        json_text = await response.text()
                        transcript_data = json.loads(json_text)
                        text = transcript_data['results']['transcripts'][0]['transcript']
                        
                        # Type the transcribed text into the active window
                        active_window = gw.getActiveWindow()
                        if active_window and text.strip():
                            pyautogui.write(text + ' ')

        finally:
            # Cleanup
            os.unlink(temp_filename)
            await self.loop.run_in_executor(
                None,
                lambda: self.s3_client.delete_object(Bucket=self.bucket_name, Key=s3_key)
            )

    async def process_audio_queue(self):
        tasks = set()
        while self.recording or not self.audio_queue.empty():
            try:
                temp_filename = self.audio_queue.get_nowait()
                task = asyncio.create_task(self.process_audio_chunk(temp_filename))
                tasks.add(task)
                task.add_done_callback(tasks.discard)
            except queue.Empty:
                await asyncio.sleep(0.1)
            
        # Wait for remaining tasks to complete
        if tasks:
            await asyncio.wait(tasks)

def main():
    app = TranscriptionApp()
    app.icon.run()

if __name__ == '__main__':
    main()
