import sys
import os
import asyncio
import datetime
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
        self.recording_thread = None
        self.processing_thread = None
        self.should_stop = threading.Event()
    
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
            self.should_stop.clear()
            self.recording_thread = threading.Thread(target=self.record_audio_chunks, daemon=True)
            self.processing_thread = threading.Thread(target=self.start_async_loop, daemon=True)
            self.recording_thread.start()
            self.processing_thread.start()
        else:
            self.should_stop.set()
            if self.recording_thread:
                self.recording_thread.join(timeout=1)
            if self.processing_thread:
                self.processing_thread.join(timeout=1)
        
    def quit_app(self):
        self.recording = False
        self.should_stop.set()
        if self.recording_thread:
            self.recording_thread.join(timeout=1)
        if self.processing_thread:
            self.processing_thread.join(timeout=1)
        self.icon.stop()

    def start_async_loop(self):
        asyncio.set_event_loop(self.loop)
        self.loop.run_until_complete(self.process_audio_queue())

    def record_audio_chunks(self):
        chunk_duration = 3  # seconds
        next_recording_ready = threading.Event()
        next_filename = None

        def start_recording():
            temp_file = tempfile.NamedTemporaryFile(suffix='.wav', delete=False)
            threading.Thread(
                target=record_audio,
                args=(chunk_duration, ),
                kwargs={'filename': temp_file.name},
                daemon=True
            ).start()
            return temp_file.name

        # Start first recording
        current_filename = start_recording()
        
        while not self.should_stop.is_set():
            # Immediately start next recording
            next_filename = start_recording()
            
            # Wait for current chunk duration
            time.sleep(chunk_duration)
            
            # Add current recording to queue and switch to next
            self.audio_queue.put(current_filename)
            current_filename = next_filename
            
        # Add final recording to queue
        if next_filename:
            time.sleep(chunk_duration)
            self.audio_queue.put(next_filename)

    async def process_audio_chunk(self, temp_filename):
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
        print(f"[{timestamp}] Processing new audio chunk: {temp_filename}")
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
            # Delete local WAV file after upload
            os.unlink(temp_filename)
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
            print(f"[{timestamp}] Uploaded to S3: {s3_key}")
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
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
            print(f"[{timestamp}] Started transcription job: {job_name}")
            
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
                        
                        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
                        print(f"[{timestamp}] Transcription completed: {text[:50]}...")

        finally:
            # Cleanup S3
            try:
                await self.loop.run_in_executor(
                    None,
                    lambda: self.s3_client.delete_object(Bucket=self.bucket_name, Key=s3_key)
                )
                timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S.%f")[:-3]
                print(f"[{timestamp}] Cleaned up S3 object: {s3_key}")
            except Exception as e:
                print(f"Error cleaning up S3 object: {e}")

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
