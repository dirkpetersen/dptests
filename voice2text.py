import sys
import os
import time
import threading
import boto3
import pystray
from PIL import Image
import pygetwindow as gw
import pyautogui
import sounddevice as sd
import soundfile as sf
import tempfile
from mic import record_audio

class TranscriptionApp:
    def __init__(self):
        self.recording = False
        self.client = boto3.client('transcribe')
        self.setup_tray()
        
    def setup_tray(self):
        # Create a simple icon (you may want to replace this with a proper icon file)
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
            threading.Thread(target=self.record_and_transcribe, daemon=True).start()
        
    def quit_app(self):
        self.recording = False
        self.icon.stop()

    def record_and_transcribe(self):
        while self.recording:
            # Create a temporary file for the audio
            with tempfile.NamedTemporaryFile(suffix='.wav', delete=False) as temp_file:
                temp_filename = temp_file.name
            
            # Record a 5-second chunk
            record_audio(5, filename=temp_filename)
            
            # Start transcription job
            job_name = f"transcription_{int(time.time())}"
            self.client.start_transcription_job(
                TranscriptionJobName=job_name,
                Media={'MediaFileUri': f"file://{temp_filename}"},
                MediaFormat='wav',
                LanguageCode='en-US'
            )
            
            # Wait for transcription to complete
            while True:
                status = self.client.get_transcription_job(TranscriptionJobName=job_name)
                if status['TranscriptionJob']['TranscriptionJobStatus'] in ['COMPLETED', 'FAILED']:
                    break
                time.sleep(0.5)
            
            # Get transcription results
            if status['TranscriptionJob']['TranscriptionJobStatus'] == 'COMPLETED':
                result = status['TranscriptionJob']['Transcript']['TranscriptFileUri']
                text = result['Results'][0]['Alternatives'][0]['Transcript']
                
                # Type the transcribed text into the active window
                active_window = gw.getActiveWindow()
                if active_window:
                    pyautogui.write(text)
            
            # Cleanup
            os.unlink(temp_filename)
            
            # Small delay before next recording
            time.sleep(0.1)

def main():
    app = TranscriptionApp()
    app.icon.run()

if __name__ == '__main__':
    main()
