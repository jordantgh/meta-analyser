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
import re

logging.basicConfig(filename="app.log", level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

load_dotenv()

Entrez.email = os.environ.get("ENTEREZ_EMAIL")
Entrez.tool = "PMC Search and Summarise"

pattern = re.compile(r'(\b(S(\d+))\s+(table|appendix|file)\b)|(\b(table|appendix|file)\s+S(\d+)\b)|(\b(supplement(ary|al)?|additional|supporting|extended|source)?\s*(info|information|file|dataset|data|table|material)\s*(supplement)?\s*(\d+)?\b)', re.IGNORECASE)


def get_supp_files(pmc_id, browser, max_retries=3):
    supp_file_dict = {}

    for attempt in range(max_retries):
        try:
            page = browser.new_page()
            page.goto(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/")

            supp_mats_tags = page.query_selector_all('.sec.suppmat')

            for sup_tag in supp_mats_tags:
                tag_id = sup_tag.get_attribute('id')
                # Get all the text outside the .sup-box
                full_text = sup_tag.inner_text()
                sup_box_text = sup_tag.query_selector('.sup-box') \
                    .inner_text() if sup_tag.query_selector('.sup-box') else ''
                
                outer_descr = full_text.replace(sup_box_text, '').strip()
                
                sup_boxes = sup_tag.query_selector_all('.sup-box')
                for sup_box in sup_boxes:
                    link = sup_box.query_selector('a')
                    href = link.get_attribute('href')
                    logging.debug(f"Link href: {href}")

                    if href.endswith('.csv') or href.endswith('.xlsx') or href.endswith('.txt'):

                        full_url = f"https://www.ncbi.nlm.nih.gov{href}"

                        # Find a useful inner description
                        inner_text = link.inner_text().strip()
                        match = pattern.search(inner_text)
                        inner_descr = match.group() if match else ""
                        
                        intext_ref = page.eval_on_selector_all(
                            f'p a[href="#{tag_id}"]',
                            '''nodes => nodes.map(node => {
                                let parent = node.closest("p");
                                return parent ? parent.innerText : "";
                            })'''
                        )

                        # Prepare the description tuple
                        desc_tuple = tuple(
                            filter(None, [outer_descr, inner_descr, intext_ref])
                        )

                        if desc_tuple:
                            supp_file_dict[full_url] = desc_tuple

            page.close()
            return supp_file_dict
        except IncompleteRead as e:
            wait_time = (2 ** attempt) + random.random()
            logging.error(
                f"IncompleteRead error: {e}. Bytes missing. Retrying in {wait_time:.2f} seconds.")
            time.sleep(wait_time)
        except Exception as e:
            wait_time = (2 ** attempt) + random.random()
            logging.error(
                f"Error encountered: {e}. Retrying in {wait_time:.2f} seconds.")
            time.sleep(wait_time)
    logging.error(f"Failed to scrape after {max_retries} attempts.")
    return supp_file_dict


def query_pmc(query, callback=None, thread=None):
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
            handle = Entrez.esearch(
                db="pmc", term=query, retstart=retstart, retmax=retmax)
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
                    handle = Entrez.efetch(
                        db="pubmed", id=pmid, rettype="medline", retmode="xml")
                    abstract_list = Entrez.read(handle, validate=False)[
                        'PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
                    abstract_text = ' '.join(
                        [str(element) for element in abstract_list])
                    handle.close()
                except (KeyError, IndexError, HTTPError):
                    abstract_text = "No abstract available."

                supp_file_dict = get_supp_files(pmcid, browser)
                supplementary_files = list(supp_file_dict.keys())
                print(supp_file_dict.values())


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
                        progress = int(
                            ((batch_num + (index + 1) / len(summaries)) / total_batches) * 100)
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
