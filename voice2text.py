import sys
import boto3
import sounddevice as sd
import numpy as np
import pystray
from PIL import Image
import pyperclip
import keyboard
from io import BytesIO
import threading
import queue
import time
from pynput.keyboard import Controller

class VoiceTranscriber:
    def __init__(self):
        self.transcribe_client = boto3.client('transcribe')
        self.keyboard = Controller()
        self.recording = False
        self.audio_queue = queue.Queue()
        self.sample_rate = 16000
        
    def start_recording(self):
        self.recording = True
        self.stream = sd.InputStream(
            channels=1,
            samplerate=self.sample_rate,
            callback=self.audio_callback
        )
        self.stream.start()

    def stop_recording(self):
        self.recording = False
        self.stream.stop()
        self.stream.close()

    def audio_callback(self, indata, frames, time, status):
        if self.recording:
            self.audio_queue.put(indata.copy())

    def process_audio(self):
        while self.recording:
            try:
                audio_chunk = self.audio_queue.get(timeout=1)
                # Convert audio to bytes
                audio_data = BytesIO()
                np.save(audio_data, audio_chunk)
                
                # Send to Amazon Transcribe
                response = self.transcribe_client.start_stream_transcription(
                    LanguageCode='en-US',
                    MediaSampleRateHertz=self.sample_rate,
                    MediaEncoding='pcm'
                )

                # Get transcription results
                for event in response['TranscriptResults']:
                    if not event['IsPartial']:
                        text = event['Alternatives'][0]['Transcript']
                        # Type the text
                        self.keyboard.type(text + ' ')

            except queue.Empty:
                continue

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
