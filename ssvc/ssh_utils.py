import os, subprocess, getpass

KEY_NAME = "ssvc"
KEY_PATH = os.path.expanduser(f"~/.ssh/{KEY_NAME}")
SSH_OPTIONS = "-o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no -o LogLevel=ERROR"

def setup_ssh_key(hostname):
    if os.path.exists(KEY_PATH):
        print("SSH key 'ssvc' already exists. Remove the existing key to repeast the setup.")
        return True

    # Generate key pair
    subprocess.run(["ssh-keygen", "-f", KEY_PATH, "-N", ""])  # No passphrase

    # Copy public key to login node 
    password = getpass.getpass(f"Password for {hostname}: ") 
    copy_cmd = ["ssh-copy-id", f"-i", KEY_PATH, f"{hostname}"]

    try:
        subprocess.run(copy_cmd, input=password.encode() + b'\n', check=True)
        print("SSH key setup successful.")
        return True
    except subprocess.CalledProcessError:
        print("SSH key setup failed.")
        return False

class SSHConnection:
    def __init__(self, hostname, username=None):
        self.hostname = hostname
        self.username = username or getpass.getuser()

    def __enter__(self):  
        self.process = subprocess.Popen(
            ["ssh", SSH_OPTIONS, self.hostname, "-l", self.username, "-i", KEY_PATH],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True  # Work with text 
        )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.process.terminate()

    def setup_port_forwarding(self, local_port, remote_host, remote_port):
        command = f" -L {local_port}:localhost:{remote_port} {remote_host}"
        print('setting up port forwarding:', command)
        self.process.stdin.write(command + "\n")
        self.process.stdin.flush() 