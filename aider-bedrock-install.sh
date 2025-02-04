#!/bin/bash

# This script installs aider.chat on Linux or OSX and sets the default LLM to Bedrock 
# and the editor to VI in multiline mode. Hit ALT+Enter instead of just Enter


# Add ~/.local/bin to PATH if not already present

RC_FILE="$HOME/.bashrc"
if [[ -n "$ZSH_VERSION" ]]; then
    RC_FILE="$HOME/.zshrc"
fi
PATH_ADDED=0
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    export PATH="$HOME/.local/bin:$PATH"
    PATH_ADDED=1
    
    if ! grep -q "export PATH=\"\$HOME/.local/bin:\$PATH\"" "$RC_FILE"; then
        echo -e "\n# Add ~/.local/bin to PATH\nexport PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$RC_FILE"
    fi
fi

# Install aider
curl -LsSf https://aider.chat/install.sh | sh

# Install boto3 using uv tool
uv tool run --from aider-chat pip install boto3

# Download and configure .aider.conf.yml if it doesn't exist
if [[ ! -f "$HOME/.aider.conf.yml" ]]; then
    curl -o "$HOME/.aider.conf.yml" https://raw.githubusercontent.com/Aider-AI/aider/refs/heads/main/aider/website/assets/sample.aider.conf.yml
    
    # Add custom settings to the top of the file
    TEMP_FILE=$(mktemp)
    cat > "$TEMP_FILE" << EOF
env-file: $HOME/.aider.env
multiline: true
dark-mode: true
cache-prompts: true
show-model-warnings: false
vim: true

EOF
    cat "$HOME/.aider.conf.yml" >> "$TEMP_FILE"
    mv "$TEMP_FILE" "$HOME/.aider.conf.yml"

# Create .aider.env file
cat > "$HOME/.aider.env" << 'EOF'
AIDER_MODEL=bedrock/us.anthropic.claude-3-5-sonnet-20241022-v2:0
# AWS_PROFILE=bedrock
# AWS_ACCESS_KEY_ID=
# AWS_SECRET_ACCESS_KEY=
# AWS_REGION=us-west-2
# AIDER_ARCHITECT=true
# AIDER_MODEL=r1
# AIDER_EDITOR_MODEL=bedrock/us.anthropic.claude-3-5-sonnet-20241022-v2:0
EOF

fi

if [[ $PATH_ADDED -eq 1 ]]; then
    echo ""	
    echo "NOTE: ~/.local/bin has been added to your PATH."
    echo "Please either:"
    echo "  1. Log out and log back in, or"
    echo "  2. Run: source $RC_FILE"
    echo "for the PATH changes to take effect."    
fi
echo "Installation complete!"
