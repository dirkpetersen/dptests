#!/bin/bash

# Add ~/.local/bin to PATH if not already present
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    export PATH="$HOME/.local/bin:$PATH"
    
    # Add to .bashrc if not already there
    if ! grep -q "export PATH=\"\$HOME/.local/bin:\$PATH\"" "$HOME/.bashrc"; then
        echo -e "\n# Add ~/.local/bin to PATH\nexport PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$HOME/.bashrc"
    fi
fi

# Install aider
curl -LsSf https://aider.chat/install.sh | sh

# Install boto3 using uv tool
uv tool run --from aider-chat pip install boto3

# Download and configure .aider.conf.yml if it doesn't exist
if [ ! -f "$HOME/.aider.conf.yml" ]; then
    curl -o "$HOME/.aider.conf.yml" https://raw.githubusercontent.com/Aider-AI/aider/refs/heads/main/aider/website/assets/sample.aider.conf.yml
    
    # Add custom settings to the top of the file
    TEMP_FILE=$(mktemp)
    cat > "$TEMP_FILE" << 'EOF'
env-file: ${HOME}/.aider.env
multiline: true
vim: true
dark-mode: true
cache-prompts: true
show-model-warnings: false

EOF
    cat "$HOME/.aider.conf.yml" >> "$TEMP_FILE"
    mv "$TEMP_FILE" "$HOME/.aider.conf.yml"
fi

# Create .aider.env file
cat > "$HOME/.aider.env" << 'EOF'
AIDER_MODEL=bedrock/us.anthropic.claude-3-5-sonnet-20241022-v2:0
# AWS_ACCESS_KEY_ID=
# AWS_SECRET_ACCESS_KEY=
# AWS_REGION=us-west-2
# AWS_PROFILE=bedrock
# AIDER_ARCHITECT=true
# AIDER_MODEL=r1
# AIDER_EDITOR_MODEL=bedrock/us.anthropic.claude-3-5-sonnet-20241022-v2:0
EOF

echo "Installation complete!"
