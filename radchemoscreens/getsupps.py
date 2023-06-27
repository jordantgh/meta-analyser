import requests
from playwright.sync_api import sync_playwright

def download_file(url):
    local_filename = url.split('/')[-1]
    # Check the file extension before downloading
    if not local_filename.endswith(('.zip', '.txt', '.csv', '.xlsx')):
        print(f"Skipping {local_filename} because it's not a .zip, .txt, .csv or .xlsx file.")
        return
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

def scrape_pmc(pmc_id):
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()

        # Navigate to the page
        page.goto(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/")

        # Find the 'data-suppmats' tag
        supp_mats_tag = page.query_selector('#data-suppmats')

        if supp_mats_tag is None:
            print(f"No supplemental materials found for {pmc_id}")
            return

        # Find all links within the 'data-suppmats' tag
        links = supp_mats_tag.query_selector_all('a')

        for link in links:
            # Extract href from each link and create full URL
            href = link.get_attribute('href')
            full_url = f"https://www.ncbi.nlm.nih.gov{href}"
            print(f"Found supplemental file link: {full_url}")

            # Download the file
            download_file(full_url)

        # Close the browser
        browser.close()

# Enter your PMC ID here
scrape_pmc("PMC7155769")