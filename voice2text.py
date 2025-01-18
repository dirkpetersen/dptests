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
import wave
import io
import uuid
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

    def transcribe(self):
        chunks = []
        while self.recording:
            try:
                chunk = self.audio_queue.get(timeout=1)
                chunks.append(chunk)
                
                # Process every 5 seconds of audio
                if len(chunks) >= (self.sample_rate * 5) // self.chunk:
                    self._process_audio_chunk(chunks)
                    chunks = []
                    
            except queue.Empty:
                continue
            except Exception as e:
                print(f"Error in transcription: {e}")
                    
    def _process_audio_chunk(self, chunks):
        # Create a temporary WAV file
        temp_filename = f"temp_{uuid.uuid4()}.wav"
        with wave.open(temp_filename, 'wb') as wf:
            wf.setnchannels(self.channels)
            wf.setsampwidth(self.p.get_sample_size(self.format))
            wf.setframerate(self.sample_rate)
            wf.writeframes(b''.join(chunks))
        
        try:
            # Start transcription job
            job_name = f"transcribe_{uuid.uuid4()}"
            self.transcribe_client.start_transcription_job(
                TranscriptionJobName=job_name,
                Media={'MediaFileUri': f"file://{os.path.abspath(temp_filename)}"},
                MediaFormat='wav',
                LanguageCode='en-US'
            )
            
            # Wait for completion
            while True:
                status = self.transcribe_client.get_transcription_job(TranscriptionJobName=job_name)
                if status['TranscriptionJob']['TranscriptionJobStatus'] in ['COMPLETED', 'FAILED']:
                    break
                time.sleep(0.5)
            
            # Get results
            if status['TranscriptionJob']['TranscriptionJobStatus'] == 'COMPLETED':
                transcript = status['TranscriptionJob']['Transcript']['Results'][0]['Alternatives'][0]['Transcript']
                self.keyboard.type(transcript + ' ')
                
        except Exception as e:
            print(f"Error processing chunk: {e}")
        finally:
            # Cleanup
            try:
                os.remove(temp_filename)
                self.transcribe_client.delete_transcription_job(TranscriptionJobName=job_name)
            except:
                pass

    def process_audio(self):
        self.transcribe()

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
