#! /usr/bin/env python3


import requests
from bs4 import BeautifulSoup
import html2text
import os

# Base URL
base_url = "https://hpc.nih.gov/apps"

def get_app_links():
    response = requests.get(base_url)
    soup = BeautifulSoup(response.content, 'lxml')  # Use lxml as the parser
    links = [a['href'] for a in soup.find_all('a', href=True) if a['href'].endswith('.html')]
    unique_links = list(set(links))  # Ensure uniqueness
    return unique_links

def fetch_and_save_app_docs():
    links = get_app_links()
    if not os.path.exists('app_docs'):
        os.makedirs('app_docs')

    for link in links:
        app_url = base_url + link
        response = requests.get(app_url)
        soup = BeautifulSoup(response.content, 'lxml')  # Use lxml as the parser
        text = soup.get_text()
        text = text.replace('sinteractive', 'grabnode')
        markdown = html2text.html2text(text)

        file_name = link.split('/')[-1]  # Extract file name from URL
        if file_name:
            file_path = os.path.join('app_docs', file_name)
            with open(file_path, 'w') as file:
                file.write(markdown)
        else:
            print(f"Warning: No valid file name extracted for URL {app_url}")

fetch_and_save_app_docs()


