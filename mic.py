import sounddevice as sd
import soundfile as sf
import argparse

def record_audio(duration, samplerate=44100, channels=2, filename='recording.wav'):
    """Record audio from default microphone"""
    # Get default input device info
    device_info = sd.query_devices(kind='input')
    print(f"Using input device: {device_info['name']}")
    print(f"Recording for {duration} seconds...")
    recording = sd.rec(
        int(duration * samplerate),
        samplerate=samplerate,
        channels=channels,
        dtype='float64'
    )
    sd.wait()  # Wait until recording is finished
    print("Recording finished")
    
    # Save to file
    sf.write(filename, recording, samplerate)
    print(f"Saved to {filename}")

def main():
    parser = argparse.ArgumentParser(description='Record audio from default microphone')
    parser.add_argument('--duration', type=float, default=5.0,
                      help='Recording duration in seconds (default: 5.0)')
    parser.add_argument('--filename', type=str, default='recording.wav',
                      help='Output filename (default: recording.wav)')
    args = parser.parse_args()

    record_audio(args.duration, filename=args.filename)

if __name__ == '__main__':
    main()
