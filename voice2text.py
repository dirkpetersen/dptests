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
        self.should_stop = threading.Event()
        
        # Audio recording settings
        self.samplerate = 44100
        self.channels = 1
        self.chunk_duration = 3  # seconds
        self.max_file_size = 100 * 1024 * 1024  # 100MB in bytes
        self.current_wav = None
        self.wav_writer = None
        self.audio_buffer = []
        self.buffer_lock = threading.Lock()
    
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
            self.audio_buffer = []  # Clear the buffer
            # Start recording thread
            threading.Thread(target=self.continuous_record, daemon=True).start()
            # Start processing thread
            threading.Thread(target=self.start_async_loop, daemon=True).start()
        else:
            self.should_stop.set()
        
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

    def create_new_wav_file(self):
        """Create a new WAV file for recording"""
        if self.current_wav:
            self.wav_writer.close()
        
        temp_file = tempfile.NamedTemporaryFile(suffix='.wav', delete=False)
        self.current_wav = temp_file.name
        self.wav_writer = sf.SoundFile(
            self.current_wav, 
            mode='w',
            samplerate=self.samplerate,
            channels=self.channels,
            format='WAV'
        )
        return self.current_wav

    def continuous_record(self):
        """Record audio continuously with a callback"""
        def audio_callback(indata, frames, time, status):
            if status:
                print(f"Status: {status}")
            
            with self.buffer_lock:
                # Write to current WAV file
                if self.wav_writer is None or self.wav_writer.tell() >= self.max_file_size:
                    self.create_new_wav_file()
                
                self.wav_writer.write(indata)
                self.audio_buffer.extend(indata[:, 0])  # Only take first channel
                
                # Check if we have enough data for a chunk
                chunk_size = int(self.samplerate * self.chunk_duration)
                while len(self.audio_buffer) >= chunk_size:
                    # Extract chunk
                    chunk = self.audio_buffer[:chunk_size]
                    self.audio_buffer = self.audio_buffer[chunk_size//2:]  # 50% overlap
                    
                    # Save chunk to temporary file
                    with tempfile.NamedTemporaryFile(suffix='.wav', delete=False) as temp_file:
                        sf.write(temp_file.name, chunk, self.samplerate)
                        self.audio_queue.put(temp_file.name)

        try:
            with sd.InputStream(
                channels=self.channels,
                samplerate=self.samplerate,
                callback=audio_callback
            ):
                self.create_new_wav_file()  # Create first WAV file
                self.should_stop.wait()
        finally:
            if self.wav_writer:
                self.wav_writer.close()

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
