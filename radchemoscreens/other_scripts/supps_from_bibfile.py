import os
import re
import requests
import bibtexparser

# Parse the BibTeX file
with open("chemo_radio_resistance_screens.bib") as bibtex_file:
    bib_database = bibtexparser.load(bibtex_file)

# Extract PMID and BibTeX key
bib_dict = {entry['pmid']: entry['ID'] for entry in bib_database.entries}

# Parse the Markdown file
with open("pubmedpmc.md") as md_file:
    md_content = md_file.read()

# Find all entries with PMID and links
entries = re.findall(r'\[(.*?)Links:\s(.*?)\]', md_content, re.DOTALL)

# Create a download directory if it does not exist
os.makedirs('downloads', exist_ok=True)

# Download the files
for entry in entries:
    pmid = re.search(r'PMID:\s(\d+)', entry[0])
    if pmid:
        pmid = pmid.group(1)
        links = re.findall(r'(https://\S+)', entry[1])
        for link in links:
            filename = link.split('/')[-1].split('?')[0]  # Original filename
            # New filename: {bibtex_key}_{PMID}_{original_filename_with ext}
            new_filename = f"{bib_dict.get(pmid, 'unknown')}_{pmid}_{filename}"
            new_filepath = os.path.join('downloads', new_filename)
            # Download the file
            r = requests.get(link, allow_redirects=True)
            open(new_filepath, 'wb').write(r.content)
            print(f"Downloaded {new_filename}")
