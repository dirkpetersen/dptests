#!/bin/bash
set -e

# Install system requirements
sudo apt-get update
sudo apt-get install -y python3-venv git python3-pip certbot traefik

# Create deployment user and structure
sudo useradd --system deployer --create-home
sudo usermod -aG sudo deployer
sudo mkdir -p /var/www
sudo chown deployer:deployer /var/www

# Configure Traefik
sudo tee /etc/traefik/traefik.yaml >/dev/null <<'EOL'
entryPoints:
  web:
    address: ":80"
  websecure:
    address: ":443"

providers:
  file:
    directory: "/etc/traefik/conf.d"
    watch: true

certificatesResolvers:
  letsencrypt:
    acme:
      email: admin@example.com
      storage: /etc/traefik/acme.json
      httpChallenge:
        entryPoint: web
EOL

sudo mkdir -p /etc/traefik/conf.d
sudo touch /etc/traefik/acme.json
sudo chmod 600 /etc/traefik/acme.json

# Configure systemd for Traefik
sudo systemctl enable traefik
sudo systemctl start traefik

# GitHub Runner setup
curl -o actions-runner-linux-x64-2.309.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.309.0/actions-runner-linux-x64-2.309.0.tar.gz
tar xzf ./actions-runner-linux-x64-2.309.0.tar.gz
rm ./actions-runner-linux-x64-2.309.0.tar.gz
./config.sh --url https://github.com/YOUR_ORG/YOUR_REPO --token YOUR_TOKEN
sudo ./svc.sh install
sudo ./svc.sh start

# Create deployment script template
sudo tee /usr/local/bin/deploy-app >/dev/null <<'EOL'
#!/bin/bash
APP_NAME=$1
APP_REPO=$2

# Clone/update repository
mkdir -p /var/www/$APP_NAME
cd /var/www/$APP_NAME
git clone $APP_REPO || git pull origin main

# Create virtual environment
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Create systemd service
sudo tee /etc/systemd/system/$APP_NAME.service >/dev/null <<EOF
[Unit]
Description=$APP_NAME Service
After=network.target

[Service]
User=deployer
WorkingDirectory=/var/www/$APP_NAME
ExecStart=/var/www/$APP_NAME/venv/bin/python -m your_app_module
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Configure Traefik route
sudo tee /etc/traefik/conf.d/$APP_NAME.yaml >/dev/null <<EOF
http:
  routers:
    ${APP_NAME}-router:
      rule: "Host(\`${APP_NAME}.s.ai.oregonstate.edu\`)"
      service: ${APP_NAME}-service
      tls:
        certResolver: letsencrypt
  
  services:
    ${APP_NAME}-service:
      loadBalancer:
        servers:
          - url: "http://localhost:8000"
EOF

sudo systemctl daemon-reload
sudo systemctl enable $APP_NAME
sudo systemctl restart $APP_NAME
EOL

sudo chmod +x /usr/local/bin/deploy-app
