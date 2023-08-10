import json
from Bio import Entrez
import time
from urllib.error import HTTPError

# Add email and script name
Entrez.email = "jtaylor81284@hotmail.com"
Entrez.tool = "PMC Search and Summarise"

def search_and_summarize_pmc(query, callback=None):
    
    # First, get the total number of results for the query
    handle = Entrez.esearch(db="pmc", term=query)
    record = Entrez.read(handle)
    total_results = int(record["Count"])
    handle.close()
    
    articles = []
    retmax = 100  # Number of results to fetch in one request
    total_batches = (total_results + retmax - 1) // retmax  # Calculate how many batches are there
    
    for batch_number, retstart in enumerate(range(0, total_results, retmax)):
        handle = Entrez.esearch(db="pmc", term=query, retstart=retstart, retmax=retmax)
        record = Entrez.read(handle)
        pmids = record["IdList"]
        handle.close()

        handle = Entrez.esummary(db="pmc", id=",".join(pmids))
        summaries = Entrez.read(handle)
        handle.close()

        for index, summary in enumerate(summaries):
            pmid = summary["ArticleIds"]["pmid"]

            try:
                handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml") 
                abstract_list = Entrez.read(handle, validate=False)['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText']
                abstract_text = ' '.join([str(element) for element in abstract_list])
                handle.close()
            except (KeyError, IndexError, HTTPError):
                abstract_text = "No abstract available."

            article_data = {
                "Title": summary["Title"],
                "Authors": summary["AuthorList"],
                "PMID": pmid,
                "PMCID": summary["ArticleIds"]["pmcid"],
                "Abstract": abstract_text
            }

            # Send the article data to the callback function (if provided)
            if callback:
                progress = int(((batch_number + (index + 1) / len(summaries)) / total_batches) * 100)
                callback(article_data, progress)
            
            # Add the article data to the articles list (this step may not be necessary anymore since we're emitting articles as they come)
            articles.append(article_data)

            time.sleep(0.5)

    return articles

def main(query_file):
    with open(query_file, 'r') as f:
        query = f.read().strip()

    articles = []
    if query:
        print("Searching and summarizing PMC with the provided query...")
        articles.extend(search_and_summarize_pmc(query))

    # Save to structured JSON
    with open("output_articles_5.json", 'w') as out:
        json.dump(articles, out, indent=4)

    print(f"Found {len(articles)} articles. Saved to output_articles.json.")

if __name__ == "__main__":
    query_file = "../pmc_searches.txt"
    main(query_file)