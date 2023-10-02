#! /bin/bash

# installs the redis compatible KV store kvrocks on Ubuntu 22.04

export DEBIAN_FRONTEND=noninteractive

sudo apt install -y git build-essential cmake libtool python3 libssl-dev
sudo mkdir -p /etc/kvrocks

mkdir -p ~/gh
cd ~/gh
git clone https://github.com/apache/kvrocks.git
cd kvrocks
./x.py build -DENABLE_OPENSSL=ON

sudo ln -s ~/gh/kvrocks/build/kvrocks /usr/local/bin/kvrocks
sudo ln -s ~/gh/kvrocks/build/kvrocks2redis /usr/local/bin/kvrocks2redis
sudo ln -s ~/gh/kvrocks/utils/systemd/kvrocks.service /etc/systemd/system/kvrocks.service
sudo cp ~/gh/kvrocks/kvrocks.conf /etc/kvrocks/kvrocks.conf.org
if ! [[ -e /etc/kvrocks/kvrocks.conf ]]; then
  sudo cp ~/gh/kvrocks/kvrocks.conf /etc/kvrocks/kvrocks.conf
fi
sudo cp ~/gh/kvrocks/utils/kvrocks2redis/kvrocks2redis.conf /etc/kvrocks/kvrocks2redis.conf.org
if ! [[ -e /etc/kvrocks/kvrocks2redis.conf ]]; then
  sudo cp ~/gh/kvrocks/utils/kvrocks2redis/kvrocks2redis.conf /etc/kvrocks/kvrocks2redis.conf
fi

sudo systemctl daemon-reload
sudo systemctl enable kvrocks.service
echo 'run: systemctl start kvrocks.service'
