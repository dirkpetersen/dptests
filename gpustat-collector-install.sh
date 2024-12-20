#! /bin/bash

mkdir -p /opt/gpustat
cd /opt/gpustat
python3 -m pip install --upgrade gpustat
curl -LO https://raw.githubusercontent.com/dirkpetersen/dptests/refs/heads/main/gpustat-collector.py
chmod +x gpustat-collector.py
crontab -l >crontab.tmp
echo "*/5 * * * * PATH=~/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin bash -c \"cd /opt/gpustat && ./gpustat-collector.py >>/dev/null\"" >>crontab.tmp
crontab crontab.tmp
rm crontab.tmp


