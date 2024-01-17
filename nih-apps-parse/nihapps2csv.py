#! /usr/bin/env python3

from bs4 import BeautifulSoup, Comment
import csv, sys, re
import requests

url_to_scrape = "https://hpc.nih.gov/apps/"  # Replace with your target URL

def parse_html_and_create_csv(html_file):
    # Read the HTML file
    with open(html_file, 'r', encoding='utf-8') as file:
        html_content = file.read()

    # Parse the HTML content
    html_content = html_content.replace(
        '<!-- End content area - do not edit below this line -->', 
        '<div class="heading"><a name="dummy" id="dummy">Dummy Category</a></div> <a href="/apps/" style="font-size:12px">back to top</a><br /> <p></p>', 1)
    soup = BeautifulSoup(html_content, 'html.parser')

    # Prepare CSV data
    csv_data = []

    # # Regex pattern to match software name and version
    # pattern = re.compile(r"(.+?)\(([^)]+)\)")

    # categories = soup.find_all('div', class_='heading')
    # for category in categories:
    #     software_category = category.get_text(strip=True)
    #     next_element = category.find_next_sibling()
    #     print(software_category, next_element)

    #     while next_element and (next_element.name != 'div' or 'heading' not in next_element.get('class', [])):
    #         if next_element.name == 'div' and 'subheading' in next_element.get('class', []):
    #             software_text = next_element.get_text(strip=True)
    #             match = pattern.search(software_text)
    #             if match:
    #                 software_name, software_version = match.groups()
    #             else:
    #                 software_name = software_text
    #                 software_version = 'Unknown'

    #             csv_data.append([software_name.lower(), software_name, software_version, software_category])

    #         next_element = next_element.find_next_sibling()

    pattern = re.compile(r"(.+?)\(([^)]+)\)")

    # Find all software categories
    categories = soup.find_all('div', class_='heading')

    for category in categories:
        software_category = category.get_text(strip=True)
        current_element = category.next_sibling

        while current_element and (current_element.name != 'div' or 'heading' not in current_element.get('class', [])):
            if current_element.name == 'div' and 'subheading' in current_element.get('class', []):
                software_text = current_element.get_text(strip=True)
                match = pattern.search(software_text)
                if match:
                    software_name, software_version = match.groups()
                else:
                    software_name = software_text
                    software_version = 'Unknown'

                csv_data.append([software_name.lower().strip(), software_name.strip(), software_version.strip(), software_category.strip()])

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
