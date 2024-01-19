#! /usr/bin/env python3
import requests
import json

class BigBadClass:
    @staticmethod
    def get_conda_package_info(channel, package_name):
        url = f'https://api.anaconda.org/package/{channel}/{package_name}'
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            required_info = {
                'channels': data.get('channels'),
                'environment_file': data.get('environment_file'),
                'extract_sources': data.get('extract_sources'),
                'install_cmd': data.get('install_cmd'),
                'install_cmds': data.get('install_cmds'),
                'prepend_to_path': data.get('prepend_to_path'),
                'remote_environment': data.get('remote_environment'),
                'requirements': data.get('requirements'),
                'staged_install': data.get('staged_install')
            }
            return required_info
        else:
            return None

# Usage
bbc = BigBadClass()
package_info = bbc.get_conda_package_info('bioconda', 'samtools')
print(json.dumps(package_info, indent=4))

