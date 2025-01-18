import sys
import boto3
import pyaudio
import numpy as np
import pystray
from PIL import Image
import win32gui
import win32con
import threading
import queue
import time
from pynput.keyboard import Controller
import websockets
import asyncio
import json

class VoiceTranscriber:
    def __init__(self):
        self.transcribe_client = boto3.client('transcribe')
        self.keyboard = Controller()
        self.recording = False
        self.audio_queue = queue.Queue()
        self.sample_rate = 16000
        self.chunk = 1024
        self.format = pyaudio.paInt16
        self.channels = 1
        self.p = pyaudio.PyAudio()
        
    def start_recording(self):
        self.recording = True
        self.stream = self.p.open(
            format=self.format,
            channels=self.channels,
            rate=self.sample_rate,
            input=True,
            frames_per_buffer=self.chunk
        )
        threading.Thread(target=self._record_audio, daemon=True).start()

    def stop_recording(self):
        self.recording = False
        if hasattr(self, 'stream'):
            self.stream.stop_stream()
            self.stream.close()

    def _record_audio(self):
        while self.recording:
            try:
                data = self.stream.read(self.chunk)
                self.audio_queue.put(data)
            except Exception as e:
                print(f"Error recording: {e}")
                break

    def process_audio(self):
        async def transcribe():
            client = boto3.client('transcribe')
            endpoint = await self._get_websocket_endpoint(client)
            
            async with websockets.connect(
                endpoint,
                ping_interval=5,
                ping_timeout=20
            ) as websocket:
                # Send configuration
                await websocket.send(json.dumps({
                    "message-type": "event",
                    "event": "start",
                    "event-type": "start_transcription",
                    "media-encoding": "pcm",
                    "sample-rate": self.sample_rate
                }))
                
                # Start audio stream thread
                audio_future = asyncio.create_task(self._stream_audio(websocket))
                
                try:
                    while self.recording:
                        msg = await websocket.recv()
                        response = json.loads(msg)
                        
                        if response.get("message-type") == "event" and response.get("event") == "transcript":
                            results = response["transcript"]["results"]
                            if results and not results[0]["is_partial"]:
                                transcript = results[0]["alternatives"][0]["transcript"]
                                self.keyboard.type(transcript + ' ')
                except Exception as e:
                    print(f"Error in transcription: {e}")
                finally:
                    audio_future.cancel()
                    
    async def _stream_audio(self, websocket):
        try:
            while self.recording:
                try:
                    chunk = self.audio_queue.get(timeout=1)
                    await websocket.send(json.dumps({
                        "message-type": "event",
                        "event": "audio-event",
                        "audio-chunk": chunk.hex()
                    }))
                except queue.Empty:
                    continue
        except Exception as e:
            print(f"Error streaming audio: {e}")
            
    async def _get_websocket_endpoint(self, client):
        response = client.start_streaming_transcription(
            LanguageCode='en-US',
            MediaEncoding='pcm',
            MediaSampleRateHertz=self.sample_rate,
            AudioStream={}
        )
        return response['WebsocketUrl']
        
        asyncio.run(transcribe())

def create_tray_icon():
    transcriber = VoiceTranscriber()
    
    def on_clicked(icon, item):
        if str(item) == 'Start':
            transcriber.start_recording()
            threading.Thread(target=transcriber.process_audio, daemon=True).start()
        elif str(item) == 'Stop':
            transcriber.stop_recording()
        elif str(item) == 'Exit':
            icon.stop()
            sys.exit()

    # Create a simple icon (you should replace this with a proper icon file)
    icon = Image.new('RGB', (64, 64), color='red')
    
    menu = pystray.Menu(
        pystray.MenuItem("Start", on_clicked),
        pystray.MenuItem("Stop", on_clicked),
        pystray.MenuItem("Exit", on_clicked)
    )
    
    tray_icon = pystray.Icon("Voice2Text", icon, "Voice2Text", menu)
    return tray_icon

def main():
    tray_icon = create_tray_icon()
    tray_icon.run()

if __name__ == "__main__":
    main()
