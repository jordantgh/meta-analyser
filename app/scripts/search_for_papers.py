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
from bs4 import BeautifulSoup

from scripts.html_sentences import (
    InsertMatcher,
    get_sentences,
    get_html_sentence
)

from unicodedata import normalize

logging.basicConfig(filename="app.log", level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

load_dotenv()

Entrez.email = os.environ.get("ENTEREZ_EMAIL")
Entrez.tool = "PMC Search and Summarise"

PATTERN = re.compile(r'(\b(S(\d+))\s+(table|appendix|file)\b)|(\b(table|appendix|file)\s+S(\d+)\b)|(\b(supplement(ary|al)?|additional|supporting|extended|source)?\s*(info|information|file|dataset|data|table|material)\s*(supplement)?\s*(\d+)?\b)', re.IGNORECASE)

BASE_URL = "https://www.ncbi.nlm.nih.gov"


def highlight_sentence_in_html(html, tag_id):
    soup = BeautifulSoup(html, 'html.parser')

    # Find the text within the <a> tags where href matches tag_id
    # e.g. <a href="#sd01">Supplementary Table 1</a>, <a href="#sd01">S1</a>
    citing_texts = []
    for a_tag in soup.find_all('a', href=f"#{tag_id}"):
        if a_tag.string:
            citing_texts.append(normalize('NFC', a_tag.string))

    if not citing_texts:
        return html  # If no citation is found, return the original html

    # Locate the parent <p> tag of the <a> tag
    paragraph = a_tag.find_parent('p') if a_tag else None
    if not paragraph:
        return html

    p_with_html = str(paragraph)
    p_textonly = paragraph.get_text()

    # Sentence tokenization and highlighting
    sentences, sentence_mappings = get_sentences(p_textonly)
    comparison = InsertMatcher(None, p_textonly, p_with_html)
    diff_ranges = comparison.get_opcodes()

    html_sentences = []
    for sentence, (start_idx, end_idx) in zip(sentences, sentence_mappings):
        html_sentence = get_html_sentence(
            sentence, start_idx, end_idx, diff_ranges, p_with_html
        )

        for cite in set(citing_texts):
            if cite in sentence:
                html_sentence = f"<strong>{html_sentence}</strong>"
                break

        html_sentences.append(html_sentence)

    paragraph.clear()
    paragraph.append(BeautifulSoup(' '.join(html_sentences), 'html.parser'))

    return str(soup)


def full_urls(description, base_url):

    def replacer(match):
        href_val = match.group(2)
        if href_val.startswith("http://") or href_val.startswith("https://"):
            # Absolute URL, no change needed
            return match.group(0)
        elif href_val.startswith("/"):
            # Relative URL, prepend domain directly
            return match.group(1) + BASE_URL + href_val
        elif href_val.startswith("#"):
            # Fragment identifier, prepend full base URL and add a slash
            return match.group(1) + base_url + href_val
        else:
            return match.group(0)

    updated_description = re.sub(
        r'(href=["\'])([^"\']+)', replacer, description)
    return updated_description


def get_supp_files(pmc_id, browser, max_retries=3):
    supp_file_dict = {}

    for attempt in range(max_retries):
        try:
            page = browser.new_page()
            page.goto(f"{BASE_URL}/pmc/articles/{pmc_id}/")

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

                        full_url = f"{BASE_URL}{href}"

                        # Find a useful inner description
                        # TODO: strip out any ':' if present
                        inner_text = link.inner_text().strip()
                        match = PATTERN.search(inner_text)
                        inner_descr = match.group() if match else ""

                        intext_ref = page.eval_on_selector_all(
                            f'p a[href="#{tag_id}"]',
                            """
                            let uniqueHTMLs = new Set();
                            nodes => nodes.map(node => {
                                let parent = node.closest("p");
                                if (parent) {
                                    let html = parent.outerHTML;
                                    if (!uniqueHTMLs.has(html)) {
                                        uniqueHTMLs.add(html);
                                        return html;
                                    }
                                }
                                return null;
                            }).filter(Boolean)
                            """
                        )

                        # Prepare the description tuple
                        outer_descr_str = f"<strong>{outer_descr}</strong>" if outer_descr else ""
                        inner_descr_str = f"<em>{inner_descr}</em>" if inner_descr else ""

                        # Formatting intext_ref entries
                        intext_ref_strs = [
                            f"<p>{highlight_sentence_in_html(ref, tag_id)}</p>" for ref in intext_ref
                        ]

                        # Joining all components for a formatted description
                        formatted_description = f"{outer_descr_str}{inner_descr_str}{''.join(intext_ref_strs)}"

                        formatted_description = full_urls(
                            formatted_description,
                            f"{BASE_URL}/pmc/articles/{pmc_id}/"
                        )

                        if formatted_description.strip():  # Checking if not just whitespace
                            supp_file_dict[full_url] = formatted_description

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
        "SupplementaryFiles": {}
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
