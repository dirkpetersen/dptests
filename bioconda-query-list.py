#! /usr/bin/env python3
import requests

class BigBadClass:
    @staticmethod
    def get_bioconda_packages():
        url = "https://api.anaconda.org/packages/bioconda"
        response = requests.get(url)
        if response.status_code == 200:
            packages = response.json()
            return [package['name'] for package in packages]
        else:
            return f"Error: {response.status_code}"

# Example Usage
packages = BigBadClass.get_bioconda_packages()
print(packages)

