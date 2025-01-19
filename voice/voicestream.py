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
        self.bytes_per_sample = 2  # 16-bit audio
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
        while not self.should_stop.is_set():
            try:
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
                    try:
                        await asyncio.gather(
                            self.receive_transcription(websocket),
                            self.send_audio(websocket)
                        )
                    except Exception as e:
                        print(f"Error in audio streaming: {e}")
                        if not self.should_stop.is_set():
                            print("Attempting to reconnect...")
                            await asyncio.sleep(1)
                            continue
                        break
            except Exception as e:
                print(f"Connection error: {e}")
                if not self.should_stop.is_set():
                    print("Attempting to reconnect...")
                    await asyncio.sleep(1)
                    continue
                break

    async def send_audio(self, websocket):
        def audio_callback(indata, frames, time, status):
            if status:
                print(f"Audio callback status: {status}")
                return  # Skip processing if there's an error
            
            # Ensure input is not empty and has correct shape
            if indata is None or indata.size == 0:
                print("Empty audio input received")
                return
            
            # Check for audio activity
            audio_level = np.abs(indata).mean()
            if audio_level > self.silence_threshold:
                self.last_audio_time = datetime.datetime.now().timestamp()
                
                # Ensure mono first
                if len(indata.shape) > 1:
                    audio_data = indata.mean(axis=1)
                else:
                    audio_data = indata.flatten()

                # Clip to [-1, 1] range
                audio_data = np.clip(audio_data, -1.0, 1.0)
                
                # Convert float32 [-1,1] to int16 [-32768,32767] with native endian
                audio_data = (audio_data * 32767).astype('int16')  # Changed from '<i2' to 'int16'
                
                # Ensure we have exactly chunk_size samples
                if len(audio_data) < self.chunk_size:
                    padding = np.zeros(self.chunk_size - len(audio_data), dtype='int16')
                    audio_data = np.concatenate([audio_data, padding])
                elif len(audio_data) > self.chunk_size:
                    audio_data = audio_data[:self.chunk_size]
                
                # Convert to raw PCM bytes (16-bit little-endian)
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
            blocksize=self.chunk_size,  # Changed from self.chunk_size * self.bytes_per_sample
            dtype=np.float32,  # Input as float32, we'll convert to int16
            latency='low'  # Reduce latency
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
                        
                    try:
                        # Create event stream message with audio data
                        headers = {
                            ':content-type': 'application/octet-stream',
                            ':event-type': 'AudioEvent',
                            ':message-type': 'event'
                        }
                        
                        # Calculate prelude and headers length
                        headers_bytes = b''
                        for key, value in headers.items():
                            header_name = key.encode('utf-8')
                            header_value = value.encode('utf-8')
                            headers_bytes += (
                                len(header_name).to_bytes(2, byteorder='big') +
                                header_name +
                                b'\x07' +  # String type
                                len(header_value).to_bytes(2, byteorder='big') +
                                header_value
                            )
                        
                        import zlib

                        # Debug headers
                        print("\nHeaders debug:")
                        for key, value in headers.items():
                            print(f"Header: {key} = {value}")
                        print(f"Headers bytes length: {len(headers_bytes)}")

                        # Create prelude (8 bytes)
                        headers_length = len(headers_bytes)
                        total_length = 8 + 4 + headers_length + len(audio_chunk) + 4  # prelude + preludeCRC + headers + payload + messageCRC
                        
                        prelude = (
                            total_length.to_bytes(4, byteorder='big') +
                            headers_length.to_bytes(4, byteorder='big')
                        )
                        
                        print("\nPrelude debug:")
                        print(f"Total message length: {total_length}")
                        print(f"Headers length: {headers_length}")
                        print(f"Prelude bytes: {list(prelude)}")

                        # Calculate prelude CRC32 (4 bytes)
                        prelude_crc = zlib.crc32(prelude) & 0xffffffff
                        prelude_with_crc = prelude + prelude_crc.to_bytes(4, byteorder='big')
                        
                        print("\nPrelude CRC debug:")
                        print(f"Prelude CRC: {prelude_crc:08x}")
                        print(f"Prelude with CRC: {prelude_with_crc.hex()}")

                        # Combine headers and payload
                        message_content = headers_bytes + audio_chunk
                        
                        print("\nPayload debug:")
                        print(f"Audio chunk length: {len(audio_chunk)}")
                        print(f"Audio chunk first 16 bytes hex: {audio_chunk[:16].hex()}")
                        print(f"Message content length: {len(message_content)}")

                        # Calculate message CRC32 (4 bytes)
                        message_crc = zlib.crc32(prelude_with_crc + message_content) & 0xffffffff
                        
                        print("\nMessage CRC debug:")
                        print(f"Message CRC: {message_crc:08x}")

                        # Combine all parts with both CRCs
                        message = prelude_with_crc + message_content + message_crc.to_bytes(4, byteorder='big')
                        
                        print("\nFinal message debug:")
                        print(f"Final message length: {len(message)}")
                        print(f"Message structure:")
                        print(f"- Prelude (8 bytes): {message[:8].hex()}")
                        print(f"- Prelude CRC (4 bytes): {message[8:12].hex()}")
                        print(f"- Headers ({headers_length} bytes): {message[12:12+headers_length].hex()}")
                        print(f"- Payload ({len(audio_chunk)} bytes): {message[12+headers_length:12+headers_length+16].hex()}...")
                        print(f"- Message CRC (4 bytes): {message[-4:].hex()}")
                        
                        print(f"Created audio event of size: {len(message)}")
                        print(f"Audio event hex: {message[:100].hex()}")
                        print(f"Audio event bytes: {list(message[:50])}")
                        
                        try:
                            await websocket.send(message)
                            print("Successfully sent audio event")
                        except websockets.exceptions.ConnectionClosed as e:
                            print(f"WebSocket connection closed while sending: {e}")
                            print(f"Close code: {e.code}")
                            print(f"Close reason: {e.reason}")
                            break
                    except Exception as e:
                        print(f"Error creating/sending audio event: {e}")
                        if "connection closed" in str(e).lower():
                            break
                        continue  # Try next chunk instead of raising
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
                    print(f"Raw response length: {len(response)}")
                    print(f"Raw response hex: {response.hex()}")
                    print(f"Raw response bytes: {list(response)}")
                    
                    try:
                        header, payload = decode_event(response)
                        print(f"Decoded header raw: {header}")
                        print(f"Decoded payload raw: {payload}")
                        print(f"Payload type: {type(payload)}")
                        
                        if isinstance(payload, bytes):
                            try:
                                payload_str = payload.decode('utf-8')
                                print(f"Decoded payload string: {payload_str}")
                                # Fix malformed JSON by adding missing braces
                                if payload_str.startswith('"Message"'):
                                    payload_str = '{' + payload_str
                                if payload_str.endswith('}'):
                                    payload = json.loads(payload_str)
                            except UnicodeDecodeError as e:
                                print(f"Unicode decode error: {e}")
                                print(f"Payload bytes: {list(payload)}")
                            except json.JSONDecodeError as e:
                                print(f"JSON decode error: {e}")
                                print(f"Attempted to parse: {payload_str}")
                        
                        print(f"Final payload: {payload}")
                        message_type = header.get(':message-type')
                        print(f"Message type: {message_type}")
                        
                        if ':exception-type' in header:
                            print(f"Exception type: {header[':exception-type']}")
                            if isinstance(payload, dict) and 'Message' in payload:
                                print(f"Exception message: {payload['Message']}")
                            
                    except Exception as e:
                        print(f"Error decoding event: {e}")
                        print(f"Exception type: {type(e)}")
                        print(f"Exception details: {str(e)}")
                        continue

                    if message_type == 'event':
                        if isinstance(payload, dict):
                            print(f"Payload keys: {payload.keys()}")
                            if 'Transcript' in payload:
                                transcript_data = payload['Transcript']
                                print(f"Transcript data: {transcript_data}")
                                results = transcript_data.get('Results', [])
                                print(f"Results: {results}")
                                if results:
                                    is_partial = results[0].get('IsPartial', True)
                                    print(f"Is partial: {is_partial}")
                                    if not is_partial:
                                        transcript = results[0]['Alternatives'][0]['Transcript']
                                        print(f"Got final transcript: {transcript}")
                                
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
