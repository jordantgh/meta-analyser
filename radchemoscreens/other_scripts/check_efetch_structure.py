from Bio import Entrez
import pprint

def inspect_record_structure_for_pmid(pmid):
    Entrez.email = "youremail@example.com"  # Always tell NCBI who you are
    
    # Use efetch to get the full details of the record with rettype set to "medline"
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
    record_details = Entrez.read(handle, validate=False)
    handle.close()
    
    return record_details

if __name__ == "__main__":
    pmid = "34624079"  # Replace with a specific pmid
    record_structure = inspect_record_structure_for_pmid(pmid)
    #pprint.pprint(record_structure)
    print(record_structure['PubmedArticle'][0]['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
