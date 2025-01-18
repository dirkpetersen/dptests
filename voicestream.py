import os
import asyncio
import websockets
import random
import string
import json
import datetime
import time
import uuid
import sounddevice as sd
import numpy as np
import pystray
from PIL import Image
import pygetwindow as gw
import pyautogui
import threading
import configparser
from pathlib import Path
from eventstream import create_audio_event, decode_event
from presigned_url import AWSTranscribePresignedURL

class TranscriptionApp:
    def load_aws_config(self):
        """Load AWS credentials and config from files or environment."""
        creds = {
            'access_key': os.getenv("AWS_ACCESS_KEY_ID", ""),
            'secret_key': os.getenv("AWS_SECRET_ACCESS_KEY", ""),
            'session_token': os.getenv("AWS_SESSION_TOKEN", ""),
            'region': os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        }
        
        # If environment variables are not set, try loading from files
        if not (creds['access_key'] and creds['secret_key']):
            config = configparser.ConfigParser()
            credentials_path = Path.home() / '.aws' / 'credentials'
            if credentials_path.exists():
                config.read(credentials_path)
                if 'default' in config:
                    creds['access_key'] = config['default'].get('aws_access_key_id', '')
                    creds['secret_key'] = config['default'].get('aws_secret_access_key', '')
                    creds['session_token'] = config['default'].get('aws_session_token', '')
        
        # Try loading region from config if not in env
        if creds['region'] == "us-east-1":  # default value
            config = configparser.ConfigParser()
            config_path = Path.home() / '.aws' / 'config'
            if config_path.exists():
                config.read(config_path)
                if 'default' in config:
                    creds['region'] = config['default'].get('region', 'us-east-1')
                elif 'profile default' in config:
                    creds['region'] = config['profile default'].get('region', 'us-east-1')
        
        if not (creds['access_key'] and creds['secret_key']):
            raise ValueError("AWS credentials not found in environment or ~/.aws/credentials")
            
        return creds
    def __init__(self):
        print("Initializing TranscriptionApp...")
        self.recording = False
        self.should_stop = asyncio.Event()
        print("Setting up system tray...")
        self.setup_tray()
        self.audio_queue = asyncio.Queue()
        
        # Audio settings
        self.sample_rate = 16000
        self.channels = 1
        self.chunk_duration = 0.1  # 100ms chunks
        self.chunk_size = int(self.sample_rate * self.chunk_duration)
        self.silence_threshold = 0.01  # Adjust this value based on testing
        self.last_audio_time = None
        
        # AWS settings
        print("Loading AWS configuration...")
        self.aws_creds = self.load_aws_config()
        print(f"Using AWS region: {self.aws_creds['region']}")
        self.transcribe_url_generator = AWSTranscribePresignedURL(
            self.aws_creds['access_key'],
            self.aws_creds['secret_key'],
            self.aws_creds.get('session_token', ''),
            self.aws_creds['region']
        )
        
        self.loop = asyncio.new_event_loop()
        
    def setup_tray(self):
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
            self.should_stop = asyncio.Event()  # Create new event
            threading.Thread(target=self.start_streaming, daemon=True).start()
        else:
            if self.loop.is_running():
                self.loop.call_soon_threadsafe(self.should_stop.set)

    def quit_app(self):
        self.recording = False
        self.should_stop.set()
        self.icon.stop()

    def start_streaming(self):
        asyncio.set_event_loop(self.loop)
        self.loop.run_until_complete(self.stream_audio())

    async def stream_audio(self):
        print("Generating presigned URL for AWS Transcribe...")
        request_url = self.transcribe_url_generator.get_request_url(
            self.sample_rate,
            "en-US",
            "pcm",
            number_of_channels=self.channels
        )
        print("Connecting to WebSocket...")

        print(f"Request URL: {request_url}")
        async with websockets.connect(request_url) as websocket:
            print("WebSocket connection established")
            await asyncio.gather(
                self.receive_transcription(websocket),
                self.send_audio(websocket)
            )

    async def send_audio(self, websocket):
        def audio_callback(indata, frames, time, status):
            if status:
                print(f"Audio callback status: {status}")
            
            # Check for audio activity
            audio_level = np.abs(indata).mean()
            if audio_level > self.silence_threshold:
                self.last_audio_time = datetime.datetime.now().timestamp()
                
                # Convert float32 to 16-bit PCM and ensure mono
                audio_data = (indata * 32767).astype(np.int16)
                if audio_data.shape[1] > 1:
                    audio_data = audio_data[:, 0]  # Take first channel if stereo
                # Ensure little-endian byte order and contiguous
                audio_data = np.ascontiguousarray(audio_data)
                audio_chunk = audio_data.tobytes()
                print(f"Audio chunk: size={len(audio_chunk)}, shape={audio_data.shape}, level={audio_level:.4f}")
                
                if len(audio_chunk) > 0:
                    try:
                        self.loop.call_soon_threadsafe(
                            self.audio_queue.put_nowait, 
                            audio_chunk
                        )
                        print(f"Queued audio chunk of size {len(audio_chunk)}")
                    except Exception as e:
                        print(f"Error queuing audio: {e}")

        stream = sd.InputStream(
            channels=self.channels,
            samplerate=self.sample_rate,
            callback=audio_callback,
            blocksize=self.chunk_size
        )
        
        with stream:
            print("Started recording...")
            while not self.should_stop.is_set():
                try:
                    audio_chunk = await asyncio.wait_for(
                        self.audio_queue.get(),
                        timeout=0.1
                    )
                    
                    # Skip sending if no recent audio activity
                    if (self.last_audio_time is None or 
                        time.time() - self.last_audio_time > 3):
                        continue
                        
                    audio_event = create_audio_event(audio_chunk)
                    await websocket.send(audio_event)
                except asyncio.TimeoutError:
                    continue
                except Exception as e:
                    print(f"Error sending audio: {e}")
                    print(f"WebSocket state: {websocket.state}")
                    break

    async def receive_transcription(self, websocket):
        last_transcript = ""
        
        while not self.should_stop.is_set():
            try:
                response = await websocket.recv()
                print(f"Received raw response length: {len(response)}")
                
                try:
                    header, payload = decode_event(response)
                    print(f"Decoded event - message type: {header.get(':message-type')}")
                    
                    if header.get(':message-type') == 'event':
                        if isinstance(payload, dict) and 'Transcript' in payload:
                            results = payload['Transcript'].get('Results', [])
                            if results and results[0].get('IsPartial', True) is False:
                                transcript = results[0]['Alternatives'][0]['Transcript']
                                print(f"Got transcript: {transcript}")
                        
                                if transcript != last_transcript:
                                    # Only type new text
                                    new_text = transcript[len(last_transcript):].strip()
                                    if new_text:
                                        print(f"Typing new text: {new_text}")
                                        active_window = gw.getActiveWindow()
                                        if active_window:
                                            pyautogui.write(new_text + ' ')
                                    last_transcript = transcript
                    elif header.get(':message-type') == 'exception':
                        error_msg = payload.get('Message', 'Unknown error')
                        print(f"Received exception: {error_msg}")
                        if "could not decode" in error_msg.lower():
                            print("Audio format error - check sample rate and encoding")
                        break
                except Exception as e:
                    print(f"Error decoding event: {e}")
                    continue
                    
            except websockets.exceptions.ConnectionClosedError:
                print("WebSocket connection closed")
                break
            except Exception as e:
                print(f"WebSocket error: {str(e)}")
                break

def main():
    app = TranscriptionApp()
    app.icon.run()

if __name__ == '__main__':
    main()
