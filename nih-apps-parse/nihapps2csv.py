#! /usr/bin/env python3

from bs4 import BeautifulSoup
import csv, sys, re
import requests

url_to_scrape = "https://hpc.nih.gov/apps/"  # Replace with your target URL

def parse_html_and_create_csv(html_file):
    # Read the HTML file
    with open(html_file, 'r', encoding='utf-8') as file:
        html_content = file.read()

    # Parse the HTML content
    soup = BeautifulSoup(html_content, 'html.parser')

    # Prepare CSV data
    csv_data = []

    # Regex pattern to match software name and version
    pattern = re.compile(r"(.+?)\(([^)]+)\)")

    # Find all software categories and iterate through them
    for category in soup.find_all('div', class_='heading'):
        software_category = category.get_text(strip=True)
        print(software_category)
        current_element = category.next_sibling

        # Traverse through siblings until the next category or end of siblings
        while current_element and (current_element.name != 'div' or 'heading' not in current_element.get('class', [])):
            if current_element.name == 'div' and 'subheading' in current_element.get('class', []):
                software_text = current_element.get_text(strip=True)
                match = pattern.search(software_text)
                if match:
                    software_name, software_version = match.groups()
                else:
                    software_name = software_text
                    software_version = 'Unknown'

                index_value = software_name.lower()
                csv_data.append([index_value.strip(), software_name.strip(), software_version.strip(), software_category.strip()])
            current_element = current_element.next_sibling

    # Write data to CSV
    with open('software_list.csv', 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['Index', 'Name', 'Version', 'Category'])
        for row in csv_data:
            writer.writerow(row)

    print(f"software_list.csv has been created successfully from {html_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        #print("Usage: python script.py <html_file>")
        response = requests.get(url_to_scrape)
        with open('index.html', 'wb') as fil:
            fil.write(response.content)
        html_file = 'index.html'
    else:
        html_file = sys.argv[1]
    parse_html_and_create_csv(html_file)
