# VS Code configs / Tips and Tricks

## Use on HPC Node via Bastion Host from Windows 

### Situation / Prerequisites 

- Secure bastion host has SSH Keys disabled and works with DUO Push 2FA
- Remote machine is WSL linux on windows with multiple entrie in ~/.ssh/config using ProxyJump to Bastion host 
- Allows multiple connections via ssh multiplexing (ControlMaster, not supported on Windows) 
- Secure bastion host has ssh key ~/.ssh/id_ed25519 for network internal communication to nodes 
- Secure bastion host has 'keychain' script in path 
- You have logged in once to bastion host since last reboot and executed this to enter the key passcode: 
  `eval $(~/.local/bin/keychain --quiet --eval id_ed25519)`

### config 

create a batch file on windows 

```
notepad %USERPROFILE%\ssh.bat
```

and paste this content 

```
%SystemRoot%\system32\wsl.exe bash -ic "eval \"\$(~/.local/bin/keychain --quiet --eval id_ed25519)\"; ssh %*"
```

open your user settings.json for VS Code

```
notepad %APPDATA%\Code\User\settings.json
```

and add 2 lines towards the end: 

```
    "remote.SSH.path": "${env:USERPROFILE}\\ssh.bat",
    "remote.SSH.useExecServer": "false",
}
```

the useExecServer is required as of Sept 2024 as there seems to be a bug in how the lastest VS code works with the remote SSH plugin 