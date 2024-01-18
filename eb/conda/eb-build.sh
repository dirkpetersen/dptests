#!/bin/bash

# Loop through each .eb file in the current directory
for file in *.eb; do
    # Check if the file exists to avoid errors in case there are no .eb files
    if [ -f "$file" ]; then
        # Run the command with the current file
        eb --rebuild "$file"
    fi
done


