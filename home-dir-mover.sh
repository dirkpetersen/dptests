#!/bin/bash

# ToDoList 
# 1) insert PATH in crontab 
# 2) echo "/arc/scratch1/dpcri/peterdir/.local/lib/python3.12/site-packages" > $HOME/.local/lib/python3.12/site-packages/custom.pth
# 3) symlink to original ~/.ssh instead of copy 
# 3) symlink to original ~/.config/systemd dir instead of copy

# Check if a destination folder is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <destination_folder>"
    exit 1
fi

# Set the destination folder
DEST=$1

# Create the destination folder if it doesn't exist
mkdir -p "$DEST"

# List of files and directories to copy
COPY_LIST=".aws .cache .config .keychain .local .ssh .vim .vscode-server bin .Xauthority .bash_history .gitconfig .lesshst .viminfo .vimrc"


# Copy files and directories
for item in $COPY_LIST; do
    echo "Copying $HOME/$item ..."
    cp -R "$HOME/$item" "$DEST/"
done

# Create symlinks for .bashrc, .bash_profile, and .profile
for file in .bashrc .bash_profile .profile; do
    if [ -f "$HOME/$file" ]; then
        mv "$HOME/$file" "$DEST/"
        ln -s "$DEST/$file" "$HOME/$file"
        echo "Created symlink for $file"
    fi
done

# Modify the .bashrc file
if [ -f "$DEST/.bashrc" ]; then
    # Use sed to replace the lines and add export HOME
    sed -i '/\. \/etc\/bashrc/{N;s@\. /etc/bashrc\nfi@. /etc/bashrc\nfi\nexport HOME='"$DEST"'@}' "$DEST/.bashrc"
    
    # Check if the modification was successful
    if grep -q "export HOME=$DEST" "$DEST/.bashrc"; then
        echo "The HOME variable has been set to $DEST in .bashrc"
    else
        echo "Warning: Could not set HOME variable in .bashrc." 
        echo "Please add this to the top of the file manually (below '. /etc/bashrc'):"
        echo "export HOME=$DEST"
    fi
else
    echo "Warning: .bashrc not found in $DEST"
fi

echo "Backup completed and symlinks created in $DEST"
echo "The HOME variable has been set to $DEST in .bashrc"
