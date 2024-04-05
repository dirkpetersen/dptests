import os
import subprocess
import getpass

KEY_NAME = "ssvc"
KEY_PATH = os.path.expanduser(f"~/.ssh/{KEY_NAME}")

def setup_ssh_key(hostname):
    if os.path.exists(KEY_PATH):
        print("SSH key 'ssvc' already exists.")
        return True

    # Generate key pair
    subprocess.run(["ssh-keygen", "-t", "rsa", "-f", KEY_PATH, "-N", ""])  # No passphrase

    # Copy public key to login node 
    password = getpass.getpass(f"Password for {hostname}: ") 
    copy_cmd = ["ssh-copy-id", f"-i", KEY_PATH, f"username@{hostname}"]

    try:
        subprocess.run(copy_cmd, input=password.encode() + b'\n', check=True)
        print("SSH key setup successful.")
        return True
    except subprocess.CalledProcessError:
        print("SSH key setup failed.")
        return False

