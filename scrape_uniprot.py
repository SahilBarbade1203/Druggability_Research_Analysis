import requests
import bs4
from bs4 import BeautifulSoup as bs
import re
import json

url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/hoxlist.txt"
page = requests.get(url)
soup = bs(page.text, 'html.parser')

text_data = soup.get_text()
file_path = "hox_proteins_data.txt"
with open(file_path, 'w', encoding='utf-8') as file:
    file.write(text_data)
# # Regular expression to capture each allergen and its corresponding protein names and IDs
# pattern = r'([A-Z][a-z]+ [a-z]+ [\d?]+[a-zA-Z]*)\s+((?:[A-Z0-9_]+\s+\([A-Z0-9]+\),?\s*)+)'
#
# # Finding all matches in the text
# matches = re.findall(pattern, text_data)
#
# # Dictionary to store allergen data
# allergen_data = {}
#
# # Iterate through the matches and populate the dictionary
# for allergen, proteins in matches:
#     allergen_data[allergen] = []
#     for protein in proteins.split(','):
#         protein_id, protein_name = re.search(r'([^\s]+)\s+\((\w+)\)', protein.strip()).groups()
#         allergen_data[allergen].append({"Protein Name": protein_name, "Protein ID": protein_id})
#
# # Appending the dictionary to a JSON file
# with open("allergen_data.json", "w") as json_file:
#     json.dump(allergen_data, json_file, indent=4)
#
# print("JSON file created successfully!")









