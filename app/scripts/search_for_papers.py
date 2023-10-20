import json
from Bio import Entrez
import time
from urllib.error import HTTPError
from playwright.sync_api import sync_playwright
import random
from dotenv import load_dotenv
import os
from http.client import IncompleteRead
import logging

logging.basicConfig(filename="app.log", level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

load_dotenv()

Entrez.email = os.environ.get("ENTEREZ_EMAIL")
Entrez.tool = "PMC Search and Summarise"

def get_supp_files(pmc_id, browser, max_retries=3):
    supplementary_files = []
    for attempt in range(max_retries):
        try:
            page = browser.new_page()
            page.goto(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/")
            supp_mats_tag = page.query_selector('#data-suppmats')
            if supp_mats_tag is not None:
                links = supp_mats_tag.query_selector_all('a')
                for link in links:
                    href = link.get_attribute('href')
                    if (href.endswith('.csv') or href.endswith('.xlsx') or href.endswith('.txt')):
                        full_url = f"https://www.ncbi.nlm.nih.gov{href}"
                        supplementary_files.append(full_url)
            page.close()
            return supplementary_files
        except IncompleteRead as e:
            wait_time = (2 ** attempt) + random.random()
            logging.error(f"IncompleteRead error: {e}. Bytes missing. Retrying in {wait_time:.2f} seconds.")
            time.sleep(wait_time)
        except Exception as e:
            wait_time = (2 ** attempt) + random.random()
            logging.error(f"Error encountered: {e}. Retrying in {wait_time:.2f} seconds.")
            time.sleep(wait_time)
    logging.error(f"Failed to scrape after {max_retries} attempts.")
    return supplementary_files

def query_pmc(query, callback=None, thread = None):
    article_json = {
        "Title": "",
        "Authors": [],
        "PMID": "",
        "PMCID": "",
        "Abstract": "",
        "SupplementaryFiles": []
    }
  
    # Launch browser once for all scrapes
    with sync_playwright() as p:
        browser = p.chromium.launch()
        handle = Entrez.esearch(db="pmc", term=query)
        record = Entrez.read(handle)
        total_results = int(record["Count"])
        handle.close()
        
        retmax = 100
        total_batches = (total_results + retmax - 1) // retmax
        
        for batch_num, retstart in enumerate(range(0, total_results, retmax)):
            if thread and thread.should_stop:
              break
            handle = Entrez.esearch(db="pmc", term=query, retstart=retstart, retmax=retmax)
            record = Entrez.read(handle)
            pmids = record["IdList"]
            handle.close()

            handle = Entrez.esummary(db="pmc", id=",".join(pmids))
            summaries = Entrez.read(handle)
            handle.close()

            for index, summary in enumerate(summaries):
                if thread and thread.should_stop:
                  break
                pmid = summary["ArticleIds"]["pmid"]
                pmcid = summary["ArticleIds"]["pmcid"]

                try:
                    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml") 
                    abstract_list = Entrez.read(handle, validate=False)['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
                    abstract_text = ' '.join([str(element) for element in abstract_list])
                    handle.close()
                except (KeyError, IndexError, HTTPError):
                    abstract_text = "No abstract available."

                supplementary_files = get_supp_files(pmcid, browser)

                article_json = {
                    "Title": summary["Title"],
                    "Authors": summary["AuthorList"],
                    "PMID": pmid,
                    "PMCID": pmcid,
                    "Abstract": abstract_text,
                    "URL": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/",
                    "SupplementaryFiles": supplementary_files
                }

                if supplementary_files:
                    if callback:
                        progress = int(((batch_num + (index + 1) / len(summaries)) / total_batches) * 100)
                        callback(article_json, progress)


                time.sleep(random.uniform(0.2, 0.5))

        browser.close()

    return article_json

def main(query_file):
    with open(query_file, 'r') as f:
        query = f.read().strip()

    articles = []
    if query:
        print("Searching and summarizing PMC with the provided query...")
        articles.extend(query_pmc(query))

    # Save to structured JSON
    with open("output_articles_6.json", 'w') as out:
        json.dump(articles, out, indent=4)

    print(f"Found {len(articles)} articles. Saved to output_articles.json.")

if __name__ == "__main__":
    query_file = "../queries/pmc_searches.txt"
    main(query_file)