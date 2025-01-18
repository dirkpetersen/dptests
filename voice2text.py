import sys
import boto3
import pyaudio
import wave
import numpy as np
import pystray
from PIL import Image
import win32gui
import win32con
import threading
import queue
import time
from pynput.keyboard import Controller
import os

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
        while self.recording:
            try:
                # Collect about 1 second of audio
                audio_data = b''
                for _ in range(int(self.sample_rate / self.chunk)):
                    chunk = self.audio_queue.get(timeout=1)
                    audio_data += chunk

                # Save temporary WAV file
                temp_wav = 'temp_audio.wav'
                with wave.open(temp_wav, 'wb') as wf:
                    wf.setnchannels(self.channels)
                    wf.setsampwidth(self.p.get_sample_size(self.format))
                    wf.setframerate(self.sample_rate)
                    wf.writeframes(audio_data)

                # Send to Amazon Transcribe
                with open(temp_wav, 'rb') as audio_file:
                    response = self.transcribe_client.start_transcription_job(
                        TranscriptionJobName=f'RealTime_{int(time.time())}',
                        Media={'MediaFileUri': f'file://{os.path.abspath(temp_wav)}'},
                        MediaFormat='wav',
                        LanguageCode='en-US'
                    )

                # Wait for transcription to complete
                job_name = response['TranscriptionJob']['TranscriptionJobName']
                while True:
                    status = self.transcribe_client.get_transcription_job(
                        TranscriptionJobName=job_name
                    )
                    if status['TranscriptionJob']['TranscriptionJobStatus'] in ['COMPLETED', 'FAILED']:
                        break
                    time.sleep(0.1)

                if status['TranscriptionJob']['TranscriptionJobStatus'] == 'COMPLETED':
                    transcript = status['TranscriptionJob']['Transcript']['Results'][0]['Alternatives'][0]['Transcript']
                    # Type the text into the active window
                    self.keyboard.type(transcript + ' ')

                # Clean up
                os.remove(temp_wav)

            except queue.Empty:
                continue
            except Exception as e:
                print(f"Error processing audio: {e}")

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
