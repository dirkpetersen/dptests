#! /usr/bin/env python3


import os
import requests
from urllib.parse import urljoin, urlparse
from bs4 import BeautifulSoup
from markdownify import markdownify as md

def download_html(url, folder_path):
    try:
        response = requests.get(url)
        if response.status_code == 200:
            # Parse the HTML content
            soup = BeautifulSoup(response.content, 'html.parser')
            
            # Find all links in the HTML content
            for link in soup.find_all('a', href=True):
                href = link.get('href')
                if href.endswith('.html'):  # Modify condition based on your need
                    full_url = urljoin(url, href)
                    download_html(full_url, folder_path)  # Recursive call for sub-pages

            # Save the HTML file
            parsed_url = urlparse(url)
            file_name = os.path.basename(parsed_url.path)
            if not file_name.endswith('.html'):
                file_name += '.html'
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'wb') as file:
                file.write(response.content)
            return file_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}")

def convert_html_to_md(html_file_path):
    with open(html_file_path, 'r') as file:
        html_content = file.read()
    markdown_content = md(html_content)
    md_file_path = html_file_path.replace('.html', '.md')
    with open(md_file_path, 'w') as file:
        file.write(markdown_content)
    print(f"Converted {html_file_path} to Markdown.")

if __name__ == "__main__":
    url_to_scrape = "https://hpc.nih.gov/apps/"  # Replace with your target URL
    local_folder = "downloaded_files"  # Local folder to save files
    if not os.path.exists(local_folder):
        os.makedirs(local_folder)

    downloaded_html_file = download_html(url_to_scrape, local_folder)
    if downloaded_html_file:
        convert_html_to_md(downloaded_html_file)


