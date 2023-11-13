#! /usr/bin/env python3

import os, sys
from markdownify import markdownify as md

def convert_html_to_md(folder_path):
    # Walk through all directories and files in the folder
    for root, dirs, files in os.walk(folder_path):
        # Filter out HTML files
        html_files = [f for f in files if f.endswith('.html')]

        # Convert each HTML file to Markdown
        for html_file in html_files:
            html_file_path = os.path.join(root, html_file)
            md_file_path = os.path.join(root, html_file.replace('.html', '.md'))

            with open(html_file_path, 'r') as file:
                html_content = file.read()

            # Convert HTML to Markdown
            markdown_content = md(html_content)

            # Write Markdown content to a new file
            with open(md_file_path, 'w') as file:
                file.write(markdown_content)

            print(f"Converted {html_file_path} to Markdown.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide a folder path.")
    else:
        folder_path = sys.argv[1]
        convert_html_to_md(folder_path)


