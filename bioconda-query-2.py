#! /usr/bin/env python3
import requests
import json

class BigBadClass:
    @staticmethod
    def get_all_conda_package_info(channel, package_name):
        url = f'https://api.anaconda.org/package/{channel}/{package_name}'
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            return None

# Usage
bbc = BigBadClass()
package_info = bbc.get_all_conda_package_info('bioconda', 'qiime2')
print(json.dumps(package_info, indent=4))

