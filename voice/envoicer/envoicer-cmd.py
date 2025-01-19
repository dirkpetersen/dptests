import os
import asyncio
import websockets
import logging
import sys
import sounddevice as sd
import numpy as np
import pygetwindow as gw
import pyautogui
import configparser
from pathlib import Path
from voice.aws_transcribe_ws.presigned_url import AWSTranscribePresignedURL
from voice.aws_transcribe_ws.eventstream import create_audio_event, decode_event

class Envoicer:
    def load_aws_config(self):
        """Load AWS credentials and config from files or environment."""
        creds = {
            'access_key': os.getenv("AWS_ACCESS_KEY_ID", ""),
            'secret_key': os.getenv("AWS_SECRET_ACCESS_KEY", ""),
            'session_token': os.getenv("AWS_SESSION_TOKEN", ""),
            'region': os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        }
        
        if not (creds['access_key'] and creds['secret_key']):
            config = configparser.ConfigParser()
            credentials_path = Path.home() / '.aws' / 'credentials'
            if credentials_path.exists():
                config.read(credentials_path)
                if 'default' in config:
                    creds['access_key'] = config['default'].get('aws_access_key_id', '')
                    creds['secret_key'] = config['default'].get('aws_secret_access_key', '')
                    creds['session_token'] = config['default'].get('aws_session_token', '')
        
        if not (creds['access_key'] and creds['secret_key']):
            raise ValueError("AWS credentials not found in environment or ~/.aws/credentials")
            
        return creds

    def __init__(self):
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Audio settings
        self.sample_rate = 16000
        self.channels = 1
        self.chunk_duration = 0.1  # 100ms chunks
        self.chunk_size = int(self.sample_rate * self.chunk_duration)
        self.bytes_per_sample = 2
        self.silence_threshold = 0.01
        
        # Setup AWS credentials
        self.aws_creds = self.load_aws_config()
        self.transcribe_url_generator = AWSTranscribePresignedURL(
            self.aws_creds['access_key'],
            self.aws_creds['secret_key'],
            self.aws_creds.get('session_token', ''),
            self.aws_creds['region']
        )
        
        self.audio_queue = asyncio.Queue()
        self.should_stop = asyncio.Event()
        self.last_transcript = ""

    async def record_audio(self):
        def audio_callback(indata, frames, time, status):
            if status:
                self.logger.warning(f"Audio status: {status}")
                return

            # Convert to mono and normalize
            audio_data = np.mean(indata, axis=1) if len(indata.shape) > 1 else indata.flatten()
            audio_data = np.clip(audio_data, -1.0, 1.0)
            
            # Convert to int16
            audio_data = (audio_data * 32767).astype(np.int16)
            
            # Queue the audio data
            audio_chunk = audio_data.tobytes()
            if len(audio_chunk) > 0:
                asyncio.run_coroutine_threadsafe(
                    self.audio_queue.put(audio_chunk), 
                    self.loop
                )

        stream = sd.InputStream(
            channels=self.channels,
            samplerate=self.sample_rate,
            callback=audio_callback,
            blocksize=self.chunk_size,
            dtype=np.float32
        )
        
        with stream:
            self.logger.info("Started recording...")
            await self.should_stop.wait()

    async def process_transcription(self, websocket):
        while not self.should_stop.is_set():
            try:
                response = await websocket.recv()
                header, payload = decode_event(response)
                
                if header[':message-type'] == 'event':
                    if isinstance(payload, dict) and 'Transcript' in payload:
                        results = payload['Transcript'].get('Results', [])
                        if results and results[0].get('Alternatives'):
                            transcript = results[0]['Alternatives'][0]['Transcript']
                            
                            # Only type new text
                            new_text = transcript[len(self.last_transcript):].strip()
                            if new_text:
                                active_window = gw.getActiveWindow()
                                if active_window:
                                    pyautogui.write(new_text + ' ')
                            self.last_transcript = transcript
                            
            except websockets.exceptions.ConnectionClosed:
                self.logger.warning("WebSocket connection closed")
                break
            except Exception as e:
                self.logger.error(f"Error processing transcription: {e}")

    async def send_audio(self, websocket):
        while not self.should_stop.is_set():
            try:
                audio_chunk = await asyncio.wait_for(self.audio_queue.get(), timeout=0.1)
                audio_event = create_audio_event(audio_chunk)
                await websocket.send(audio_event)
            except asyncio.TimeoutError:
                continue
            except Exception as e:
                self.logger.error(f"Error sending audio: {e}")
                break

    async def run(self):
        self.loop = asyncio.get_running_loop()
        
        while not self.should_stop.is_set():
            try:
                request_url = self.transcribe_url_generator.get_request_url(
                    self.sample_rate,
                    "en-US",
                    "pcm",
                    number_of_channels=self.channels
                )
                
                async with websockets.connect(request_url) as websocket:
                    self.logger.info("Connected to AWS Transcribe")
                    await asyncio.gather(
                        self.record_audio(),
                        self.send_audio(websocket),
                        self.process_transcription(websocket)
                    )
            except Exception as e:
                self.logger.error(f"Connection error: {e}")
                if not self.should_stop.is_set():
                    await asyncio.sleep(1)
                    continue
                break

    def start(self):
        try:
            asyncio.run(self.run())
        except KeyboardInterrupt:
            self.logger.info("Stopping transcription...")
            self.should_stop.set()

def main():
    envoicer = Envoicer()
    envoicer.start()

if __name__ == "__main__":
    main()
