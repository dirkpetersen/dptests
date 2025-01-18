import os
import asyncio
import websockets
import json
import datetime
import uuid
import sounddevice as sd
import numpy as np
import pystray
from PIL import Image
import pygetwindow as gw
import pyautogui
import threading
from eventstream import create_audio_event, decode_event
from presigned_url import AWSTranscribePresignedURL

class TranscriptionApp:
    def __init__(self):
        self.recording = False
        self.should_stop = threading.Event()
        self.setup_tray()
        
        # Audio settings
        self.sample_rate = 16000
        self.channels = 1
        self.chunk_duration = 0.1  # 100ms chunks
        self.chunk_size = int(self.sample_rate * self.chunk_duration)
        
        # AWS settings
        self.region = os.getenv("AWS_DEFAULT_REGION", "us-east-1")
        self.access_key = os.getenv("AWS_ACCESS_KEY_ID", "")
        self.secret_key = os.getenv("AWS_SECRET_ACCESS_KEY", "")
        self.session_token = os.getenv("AWS_SESSION_TOKEN", "")
        
        self.transcribe_url_generator = AWSTranscribePresignedURL(
            self.access_key, 
            self.secret_key,
            self.session_token,
            self.region
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
            self.should_stop.clear()
            threading.Thread(target=self.start_streaming, daemon=True).start()
        else:
            self.should_stop.set()

    def quit_app(self):
        self.recording = False
        self.should_stop.set()
        self.icon.stop()

    def start_streaming(self):
        asyncio.set_event_loop(self.loop)
        self.loop.run_until_complete(self.stream_audio())

    async def stream_audio(self):
        websocket_key = str(uuid.uuid4())
        extra_headers = {
            "Origin": "https://localhost",
            "Sec-Websocket-Key": websocket_key,
            "Sec-Websocket-Version": "13",
            "Connection": "keep-alive"
        }
        
        request_url = self.transcribe_url_generator.get_request_url(
            self.sample_rate,
            "en-US",
            "pcm",
            number_of_channels=self.channels
        )

        async with websockets.connect(
            request_url,
            extra_headers=extra_headers,
            ping_timeout=None
        ) as websocket:
            await asyncio.gather(
                self.receive_transcription(websocket),
                self.send_audio(websocket)
            )

    async def send_audio(self, websocket):
        def audio_callback(indata, frames, time, status):
            if status:
                print(f"Status: {status}")
            
            audio_chunk = indata.tobytes()
            if len(audio_chunk) > 0:
                audio_event = create_audio_event(audio_chunk)
                asyncio.run_coroutine_threadsafe(
                    websocket.send(audio_event),
                    self.loop
                )

        with sd.InputStream(
            channels=self.channels,
            samplerate=self.sample_rate,
            callback=audio_callback,
            blocksize=self.chunk_size
        ):
            await self.should_stop.wait()

    async def receive_transcription(self, websocket):
        last_transcript = ""
        
        while not self.should_stop.is_set():
            try:
                response = await websocket.recv()
                header, payload = decode_event(response)
                
                if header[':message-type'] == 'event':
                    results = payload.get('Transcript', {}).get('Results', [])
                    if results and len(results) > 0:
                        transcript = results[0]['Alternatives'][0]['Transcript']
                        
                        if transcript != last_transcript:
                            # Only type new text
                            new_text = transcript[len(last_transcript):].strip()
                            if new_text:
                                active_window = gw.getActiveWindow()
                                if active_window:
                                    pyautogui.write(new_text + ' ')
                            last_transcript = transcript
                            
                elif header[":message-type"] == 'exception':
                    print(f"Error: {payload.get('Message', 'Unknown error')}")
                    
            except websockets.exceptions.ConnectionClosedError:
                print("Connection closed")
                break
            except Exception as e:
                print(f"Error: {str(e)}")
                break

def main():
    app = TranscriptionApp()
    app.icon.run()

if __name__ == '__main__':
    main()
