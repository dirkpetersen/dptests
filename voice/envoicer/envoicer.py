import os
import asyncio
import websockets
import logging
import sys
import string
import random
import pyaudio
import wave
import keyboard
import pyautogui
import pygetwindow as gw
from presigned_url import AWSTranscribePresignedURL
from eventstream import create_audio_event, decode_event

class Envoicer:
    def __init__(self):
        # Configure logging
        logging.basicConfig(
            stream=sys.stderr,
            level=logging.DEBUG,  # Changed to DEBUG level
            format='%(asctime)s.%(msecs)03d %(levelname)s: %(message)s',
            datefmt='%H:%M:%S'
        )
        
        # Audio settings optimized for lowest latency
        self.CHANNELS = 1
        self.RATE = 16000  # Higher quality for speech recognition
        self.CHUNK = 1024  # Increased for better stability
        self.FORMAT = pyaudio.paInt16
        self.running = False
        self.last_text = ""
        self.sent_sentences = set()  # Track sent sentences
        self.partial_stability_counter = 0
        self.silence_threshold = 300  # Lower threshold for audio activity
        self.debug_audio = True  # Enable audio level debugging
        
        # AWS Configuration
        self.access_key = os.getenv("AWS_ACCESS_KEY_ID", "")
        self.secret_key = os.getenv("AWS_SECRET_ACCESS_KEY", "")
        self.session_token = os.getenv("AWS_SESSION_TOKEN", "")
        self.region = os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        
        # Initialize PyAudio
        self.audio = pyaudio.PyAudio()
        
        # Configure pyautogui
        pyautogui.FAILSAFE = True
        pyautogui.PAUSE = 0.001  # Minimal delay between actions
        
    def get_default_input_device_info(self):
        """Get and log information about the default input device"""
        try:
            default_input = self.audio.get_default_input_device_info()
            logging.info(f"Using input device: {default_input['name']}")
            logging.info(f"Device info: {default_input}")
            return default_input
        except Exception as e:
            logging.error(f"Error getting input device info: {e}")
            return None

    def is_audio_active(self, audio_data):
        """Check if there's significant audio activity"""
        audio_level = max(abs(int.from_bytes(audio_data[i:i+2], 'little', signed=True)) 
                         for i in range(0, len(audio_data), 2))
        
        if self.debug_audio and audio_level > 100:  # Only log significant levels
            logging.debug(f"Audio level: {audio_level}")
        
        return audio_level > self.silence_threshold

    async def record_and_stream(self, websocket):
        self.get_default_input_device_info()
        stream = self.audio.open(
            format=self.FORMAT,
            channels=self.CHANNELS,
            rate=self.RATE,
            input=True,
            frames_per_buffer=self.CHUNK
        )
        logging.info(f"Started recording with: {self.RATE}Hz, {self.CHANNELS} channels, chunk size: {self.CHUNK}")
        
        try:
            consecutive_errors = 0
            while self.running:
                try:
                    data = stream.read(self.CHUNK, exception_on_overflow=False)
                    if len(data) > 0:
                        if self.is_audio_active(data):
                            logging.debug("Sending audio data to AWS")
                            audio_event = create_audio_event(data)
                            await websocket.send(audio_event)
                            consecutive_errors = 0
                        await asyncio.sleep(0.001)  # Small sleep to prevent CPU overload
                except websockets.exceptions.ConnectionClosedError:
                    consecutive_errors += 1
                    if consecutive_errors > 5:
                        logging.error("Too many consecutive connection errors")
                        break
                    continue
        finally:
            stream.stop_stream()
            stream.close()

    async def receive_transcription(self, websocket):
        try:
            while self.running:
                try:
                    response = await websocket.recv()
                    header, payload = decode_event(response)
                    logging.debug(f"Received response - header: {header}")
                    logging.debug(f"Payload: {payload}")
                
                    if header[':message-type'] == 'event':
                        if 'Transcript' in payload and len(payload['Transcript']['Results']) > 0:
                            transcript = payload['Transcript']['Results'][0]
                            if 'Alternatives' in transcript and len(transcript['Alternatives']) > 0:
                                text = transcript['Alternatives'][0]['Transcript'].strip()
                                is_partial = transcript.get('IsPartial', True)
                                
                                logging.info(f"Transcribed{'(partial)' if is_partial else ''}: {text}")
                            
                            if text:
                                try:
                                    if is_partial:
                                        # Only update text if it ends with a complete word
                                        # (no word fragments at the end)
                                        if text != self.last_text and (text.endswith(' ') or text.endswith('.')):
                                            self.partial_stability_counter = 0
                                            # Delete previous text and type new text
                                            if self.last_text:
                                                for _ in range(len(self.last_text)):
                                                    pyautogui.press('backspace')
                                                logging.info(f"SEND: [BACKSPACE x{len(self.last_text)}]")
                                            pyautogui.write(text)
                                            logging.info(f"SEND: {text}")
                                            self.last_text = text
                                        else:
                                            self.partial_stability_counter += 1
                                    else:
                                        # Only send if we haven't sent this sentence before
                                        if text not in self.sent_sentences:
                                            # For final text, add a space after if it doesn't end with punctuation
                                            if self.last_text:
                                                for _ in range(len(self.last_text)):
                                                    pyautogui.press('backspace')
                                                logging.info(f"SEND: [BACKSPACE x{len(self.last_text)}]")
                                            ending = ' ' if not text.endswith(('.', '!', '?')) else ''
                                            pyautogui.write(text + ending)
                                            logging.info(f"SEND: {text}{ending}")
                                            self.sent_sentences.add(text)  # Add to sent sentences
                                            self.last_text = ""
                                            self.partial_stability_counter = 0
                                except Exception as e:
                                    logging.error(f"Error sending keys: {e}")
                except websockets.exceptions.ConnectionClosedOK:
                    logging.info("Streaming completed successfully - reconnecting...")
                    return  # Allow reconnection in connect_to_websocket
                except websockets.exceptions.ConnectionClosedError:
                    logging.error("WebSocket connection closed unexpectedly")
                    return  # Allow reconnection in connect_to_websocket
                except Exception as e:
                    logging.exception("Error in receive_transcription")
                    return  # Allow reconnection in connect_to_websocket
        except Exception as e:
            logging.exception("Fatal error in receive_transcription")

    async def connect_to_websocket(self):
        max_retries = 3
        retry_delay = 2
        attempt = 0
        
        while self.running and attempt < max_retries:
            try:
                attempt += 1
                logging.info(f"Connecting to AWS Transcribe (attempt {attempt}/{max_retries})")
                
                transcribe_url_generator = AWSTranscribePresignedURL(
                    self.access_key, self.secret_key, self.session_token, self.region
                )
                
                websocket_key = ''.join(random.choices(
                    string.ascii_uppercase + string.ascii_lowercase + string.digits, k=20
                ))
                
                headers = {
                    "Origin": "https://localhost",
                    "Sec-Websocket-Key": websocket_key,
                    "Sec-Websocket-Version": "13",
                    "Connection": "keep-alive"
                }
                
                request_url = transcribe_url_generator.get_request_url(
                    sample_rate=self.RATE,
                    language_code="en-US",
                    media_encoding="pcm",
                    number_of_channels=self.CHANNELS,
                    enable_channel_identification=False,
                    enable_partial_results_stabilization=True,
                    partial_results_stability="medium"
                )
                
                async with websockets.connect(
                    request_url,
                    additional_headers=headers,
                    ping_timeout=20,  # Add ping timeout
                    ping_interval=15,  # Keep connection alive with pings
                    close_timeout=5,
                    max_size=2**24,  # Increase max message size
                    compression=None  # Disable compression for better performance
                ) as websocket:
                    logging.info("Connected to AWS Transcribe")
                    await asyncio.gather(
                        self.record_and_stream(websocket),
                        self.receive_transcription(websocket)
                    )
                    
            except websockets.exceptions.ConnectionClosedError as e:
                if attempt < max_retries:
                    logging.warning(f"Connection closed, retrying in {retry_delay} seconds... ({e})")
                    await asyncio.sleep(retry_delay)
                else:
                    logging.error(f"Failed to maintain connection after {max_retries} attempts")
                    break
            except Exception as e:
                logging.exception(f"Unexpected error in connection: {e}")
                if attempt < max_retries:
                    await asyncio.sleep(retry_delay)
                else:
                    break

    def start(self):
        """Start the voice transcription service"""
        self.running = True
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        try:
            # Register hotkey to stop the service
            keyboard.add_hotkey('ctrl+shift+x', self.stop)
            
            print("\nEnvoicer started. Press Ctrl+Shift+X to stop.")
            print("Listening for audio input...")
            loop.run_until_complete(self.connect_to_websocket())
        except KeyboardInterrupt:
            self.stop()
        finally:
            self.audio.terminate()
            loop.close()
    
    def stop(self):
        """Stop the voice transcription service"""
        self.running = False
        print("\nEnvoicer stopped.")

def main():
    envoicer = Envoicer()
    envoicer.start()

if __name__ == "__main__":
    main()
