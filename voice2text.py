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
import asyncio
from amazon_transcribe.client import TranscribeStreamingClient
from amazon_transcribe.handlers import TranscriptResultStreamHandler
from amazon_transcribe.model import TranscriptEvent

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

    async def process_stream(self, stream):
        async client = TranscribeStreamingClient(region="us-west-2")
        
        class MyEventHandler(TranscriptResultStreamHandler):
            def __init__(self, keyboard):
                super().__init__()
                self.keyboard = keyboard
                
            async def handle_transcript_event(self, transcript_event: TranscriptEvent):
                results = transcript_event.transcript.results
                for result in results:
                    if not result.is_partial:
                        transcript = result.alternatives[0].transcript
                        self.keyboard.type(transcript + ' ')

        handler = MyEventHandler(self.keyboard)
        
        async with client.start_stream_transcription(
            language_code="en-US",
            media_sample_rate_hz=self.sample_rate,
            media_encoding="pcm"
        ) as stream:
            handler.handle_events(stream.output_stream)
            while self.recording:
                try:
                    chunk = self.audio_queue.get(timeout=1)
                    await stream.input_stream.send_audio_event(audio_chunk=chunk)
                except queue.Empty:
                    continue
                except Exception as e:
                    print(f"Error processing audio: {e}")
                    break
            await stream.input_stream.end_stream()

    def process_audio(self):
        asyncio.run(self.process_stream())

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
