import os
import asyncio
import websockets
import logging
import sys
import string
import random
import pyaudio
import wave
import win32com.client
import keyboard
from presigned_url import AWSTranscribePresignedURL
from eventstream import create_audio_event, decode_event

class Envoicer:
    def __init__(self):
        # Configure logging
        logging.basicConfig(stream=sys.stderr, level=logging.INFO)
        
        # Audio settings optimized for lowest latency
        self.CHANNELS = 1
        self.RATE = 16000  # Higher quality for speech recognition
        self.CHUNK = 512   # Smaller chunks for even lower latency
        self.FORMAT = pyaudio.paInt16
        self.running = False
        self.last_text = ""
        
        # AWS Configuration
        self.access_key = os.getenv("AWS_ACCESS_KEY_ID", "")
        self.secret_key = os.getenv("AWS_SECRET_ACCESS_KEY", "")
        self.session_token = os.getenv("AWS_SESSION_TOKEN", "")
        self.region = os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        
        # Initialize PyAudio
        self.audio = pyaudio.PyAudio()
        
        # Initialize Windows shell for typing
        self.shell = win32com.client.Dispatch("WScript.Shell")
        
    async def record_and_stream(self, websocket):
        stream = self.audio.open(
            format=self.FORMAT,
            channels=self.CHANNELS,
            rate=self.RATE,
            input=True,
            frames_per_buffer=self.CHUNK
        )
        
        try:
            while self.running:
                data = stream.read(self.CHUNK, exception_on_overflow=False)
                if len(data) > 0:
                    audio_event = create_audio_event(data)
                    await websocket.send(audio_event)
                # No sleep needed - we want to stream as fast as possible
        finally:
            stream.stop_stream()
            stream.close()

    async def receive_transcription(self, websocket):
        try:
            while self.running:
                response = await websocket.recv()
                header, payload = decode_event(response)
                
                if header[':message-type'] == 'event' and len(payload['Transcript']['Results']) > 0:
                    transcript = payload['Transcript']['Results'][0]
                    text = transcript['Alternatives'][0]['Transcript'].strip()
                    
                    if text and text != self.last_text:
                        # For partial results, backspace the previous text
                        if transcript.get('IsPartial', True) and self.last_text:
                            self.shell.SendKeys("{BS " + str(len(self.last_text)) + "}")
                        
                        # Send the new text
                        self.shell.SendKeys(text)
                        self.last_text = text
                        
                        # Add space only for final results
                        if not transcript.get('IsPartial', True):
                            self.shell.SendKeys(" ")
                            self.last_text = ""
        except websockets.exceptions.ConnectionClosedError:
            logging.error("WebSocket connection closed")
        except Exception as e:
            logging.exception("Error in receive_transcription")

    async def connect_to_websocket(self):
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
            partial_results_stability="high"
        )
        
        async with websockets.connect(
            request_url,
            additional_headers=headers,
            ping_timeout=None
        ) as websocket:
            await asyncio.gather(
                self.record_and_stream(websocket),
                self.receive_transcription(websocket)
            )

    def start(self):
        """Start the voice transcription service"""
        self.running = True
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        try:
            # Register hotkey to stop the service
            keyboard.add_hotkey('ctrl+shift+x', self.stop)
            
            print("Envoicer started. Press Ctrl+Shift+X to stop.")
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
