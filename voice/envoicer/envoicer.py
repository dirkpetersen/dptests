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
            level=logging.INFO,  # Change to INFO level to reduce debug messages
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
        self.debug_audio = False  # Disable audio level debugging
        
        # AWS Configuration
        self.access_key = os.getenv("AWS_ACCESS_KEY_ID", "")
        self.secret_key = os.getenv("AWS_SECRET_ACCESS_KEY", "")
        self.session_token = os.getenv("AWS_SESSION_TOKEN", "")
        self.region = os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        
        # Initialize PyAudio
        self.audio = pyaudio.PyAudio()
        
        # Configure pyautogui
        pyautogui.FAILSAFE = True
        pyautogui.PAUSE = 0.01  # Slightly longer delay for reliability
        self.last_active_window = None
        
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
                
                    if header[':message-type'] == 'event':
                        if 'Transcript' in payload and len(payload['Transcript']['Results']) > 0:
                            transcript = payload['Transcript']['Results'][0]
                            if 'Alternatives' in transcript and len(transcript['Alternatives']) > 0:
                                text = transcript['Alternatives'][0]['Transcript'].strip()
                                is_partial = transcript.get('IsPartial', True)
                                
                                # Get current active window
                                active_window = gw.getActiveWindow()
                                
                                # Only log window info if it changed
                                if active_window:
                                    current_window_info = f"{active_window.title} (Handle: {active_window._hWnd})"
                                    if not hasattr(self, '_last_window_info') or self._last_window_info != current_window_info:
                                        logging.info(f"Active window: {current_window_info}")
                                        self._last_window_info = current_window_info
                                elif not hasattr(self, '_last_window_info') or self._last_window_info is not None:
                                    logging.info("No active window")
                                    self._last_window_info = None
                                
                                print(f"\nTranscript: {text}")
                            
                            if text:
                                try:
                                    if is_partial:
                                        # Only update text if it ends with a complete word
                                        if text != self.last_text and (text.endswith(' ') or text.endswith('.')):
                                            self.partial_stability_counter = 0
                                            try:
                                                active_window = gw.getActiveWindow()
                                                if active_window:
                                                    # Log window details before activation
                                                    logging.info(f"Found active window: {active_window.title}")
                                                    logging.info(f"Window details - PID: {active_window._hWnd}, Size: {active_window.size}, Position: {active_window.topleft}")
                                                    
                                                    # Ensure window is active and ready
                                                    active_window.activate()
                                                    active_window.restore()
                                                    pyautogui.sleep(0.1)
                                                    
                                                    # Verify window activation
                                                    current_window = gw.getActiveWindow()
                                                    if current_window and current_window._hWnd == active_window._hWnd:
                                                        logging.info(f"Successfully activated window: {current_window.title}")
                                                    else:
                                                        logging.warning(f"Window activation may have failed - Current: {current_window.title if current_window else 'None'}")
                                                    
                                                    # Type the new text
                                                    new_text = text[len(self.last_text):].strip()
                                                    if new_text:
                                                        pyautogui.write(new_text + ' ')
                                                        logging.info(f"Typed: {new_text}")
                                                else:
                                                    # Log all windows to help debug
                                                    all_windows = gw.getAllWindows()
                                                    logging.warning("No active window found. Available windows:")
                                                    for window in all_windows:
                                                        logging.warning(f"- {window.title} (pid: {window._hWnd})")
                                            except Exception as e:
                                                logging.error(f"Error getting window info: {e}")
                                                self.last_text = text
                                        else:
                                            self.partial_stability_counter += 1
                                    else:
                                        # Only send if we haven't sent this sentence before
                                        if text not in self.sent_sentences:
                                            # For final text, add a space after if it doesn't end with punctuation
                                            active_window = gw.getActiveWindow()
                                            if active_window:
                                                # Log window details before proceeding
                                                logging.info(f"Current active window: {active_window.title}")
                                                logging.info(f"Window details - Handle: {active_window._hWnd}, Size: {active_window.size}, Position: {active_window.topleft}")
                                                    
                                                # Only proceed if window changed or no previous window
                                                if active_window != self.last_active_window:
                                                    logging.info(f"Window changed from: {self.last_active_window.title if self.last_active_window else 'None'}")
                                                    self.last_active_window = active_window
                                                    pyautogui.sleep(0.1)  # Small pause when switching windows
                                                
                                                try:
                                                    # Ensure window is active
                                                    active_window.activate()
                                                    if self.last_text:
                                                        pyautogui.hotkey('ctrl', 'backspace')  # Delete last word
                                                        logging.info(f"SEND to '{active_window.title}': [DELETE WORD]")
                                                    ending = ' ' if not text.endswith(('.', '!', '?')) else ''
                                                    pyautogui.write(text + ending, interval=0.01)  # Slightly slower typing
                                                    logging.info(f"SEND to '{active_window.title}': {text}{ending}")
                                                except Exception as e:
                                                    logging.error(f"Failed to send text: {e}")
                                            else:
                                                logging.warning("No active window found - text not sent")
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
                    ping_timeout=20,
                    ping_interval=15,
                    close_timeout=5,
                    max_size=2**24,
                    compression=None
                ) as websocket:
                    logging.info("Connected to AWS Transcribe")
                    try:
                        await asyncio.gather(
                            self.record_and_stream(websocket),
                            self.receive_transcription(websocket)
                        )
                    except websockets.exceptions.ConnectionClosedOK:
                        logging.info("Connection closed normally, reconnecting...")
                        await asyncio.sleep(1)  # Brief pause before reconnecting
                        continue
                    except websockets.exceptions.ConnectionClosedError as e:
                        if attempt < max_retries:
                            logging.warning(f"Connection closed unexpectedly, retrying in {retry_delay} seconds... ({e})")
                            await asyncio.sleep(retry_delay)
                            continue
                        else:
                            logging.error(f"Failed to maintain connection after {max_retries} attempts")
                            break
                    
            except Exception as e:
                logging.exception(f"Unexpected error in connection: {e}")
                if attempt < max_retries:
                    await asyncio.sleep(retry_delay * attempt)  # Exponential backoff
                    continue
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
