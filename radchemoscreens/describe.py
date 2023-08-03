import time
import os
import pandas as pd
import openai
import bibtexparser
import json

openai.api_key = "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"

def analyze_dataset(dataset_filepath, paper_title):
    # Load the dataset
    file_ext = os.path.splitext(dataset_filepath)[-1]
    if file_ext == '.csv':
        df = pd.read_csv(dataset_filepath)
    elif file_ext == '.txt':
        df = pd.read_csv(dataset_filepath, delimiter='\t')
    elif file_ext == '.xlsx':
        # Read all sheets
        dfs = pd.read_excel(dataset_filepath, sheet_name=None)
    else:
        print(f"Unsupported file format: {file_ext}")
        return
    
    # filname without base directory or extension
    filename = os.path.splitext(os.path.basename(dataset_filepath))[0]

    process_dataset(df, filename, paper_title)

def convert_to_markdown(json_obj):
    markdown_str = "# ---\n"
    for key, value in json_obj.items():
        if key == "sample":
            markdown_str += f"{key}:\n{value}\n\n"
        else:
            markdown_str += f"{key}: {value}\n\n"
    return markdown_str

prompt_examples = [
    {
        "role": "user",
        "content": """
        **Article:** The genomic landscape of sensitivity to arsenic trioxide uncovered by genome-wide CRISPR-Cas9 screening.
        **sample:**
        | id       |  pos|score |  pos|p-value |  pos|fdr |  pos|rank |  pos|goodsgrna |  pos|lfc |
        |:---------|------------:|--------------:|----------:|-----------:|----------------:|----------:|
        | KEAP1    |    2.87e-10 |      2.59e-07 |  0.002475 |          1 |               3 |   10.566  |
        | C19orf43 |    9.42e-07 |      2.85e-06 |  0.018152 |          2 |               3 |    8.7163 |
        | PIGA     |    2.12e-05 |      7.56e-05 |  0.294554 |          3 |               2 |    4.5552 |
        | BBC3     |    2.24e-05 |      8e-05    |  0.294554 |          4 |               4 |    7.5085 |
        | PRKCSH   |    2.64e-05 |      9.25e-05 |  0.294554 |          5 |               2 |    1.8063 |
        **comment:** Hint: this is what MAGeCK statistics look like. Generally in the form "[pos/neg]|[summary statistic]" e.g. pos|p-value, neg|goodsgrna,  pos|lfc, etc.
        """
    },
    {
        "role": "assistant",
        "content": """
        {
           "article_title": "The genomic landscape of sensitivity to arsenic trioxide uncovered by genome-wide CRISPR-Cas9 screening.",
           "dataset_filename": "chenGenomicLandscapeSensitivity2023_37251921_Table_1_ATO_normalized-total.gene_summa_t1",
           "header_row": "0",
           "sgRNA_sequence": "None",
           "gene_identifier": "HGNC gene symbol",
           "dataset_type": "Data",
           "statistic_aggregation": "Per gene",
           "metrics": ["MAGeCK statistics"],
           "conditions": ["ATO"]
        }
        """
    },
    {
        "role": "user",
        "content": """
        **Article:** Genome scale CRISPR Cas9a knockout screen reveals genes that control glioblastoma susceptibility to the alkylating agent temozolomide
        **sample:**
        | A1BG_CATCTTCTTTCACCTGAACG   |   677 |
        |:----------------------------|------:|
        | A1BG_CTCCGGGGAGAACTCCGGCG   |  1248 |
        | A1BG_TCTCCATGGTGCATCAGCAC   |    88 |
        | A1BG_TGGAAGTCCACTCCACTCAG   |   842 |
        | A1CF_ACAGGAAGAATTCAGTTATG   |   972 |
        | A1CF_AGTTATGTTAGGTATACCCG   |   749 |
        """
    },
    {
        "role": "assistant",
        "content": """
        {
          "article_title": "Genome scale CRISPR Cas9a knockout screen reveals genes that control glioblastoma susceptibility to the alkylating agent temozolomide",
          "dataset_filename": "awahGenomeScaleCRISPR2022_35990011_5953_DMSO_Day14_R2_sheet_t1.csv",
          "header_row": "None",
          "sgRNA_sequence": "Concatenated",
          "gene_identifier": "HGNC gene symbol",
          "dataset_type": "Data",
          "statistic_aggregation": "Per sgRNA",
          "metrics": ["Raw count"],
          "conditions": ["DMSO"]
        }
        """
    },
    {
        "role": "user",
        "content": """
        **Article:** Cooperation of the ATM and Fanconi Anemia/BRCA Pathways in Double-Strand Break End Resection
        **sample:**
        | Barcode Sequence     | Annotated Gene Symbol   |   Annotated Gene ID |
        |:---------------------|:------------------------|--------------------:|
        | AAAAAAAATCCGGACAATGG | SLC25A24                |               29957 |
        | AAAAAAAGGATGGTGATCAA | FASTKD3                 |               79072 |
        | AAAAAAATGACATTACTGCA | BCAS2                   |               10286 |
        | AAAAAAATGTCAGTCGAGTG | GPR18                   |                2841 |
        | AAAAAACACAAGCAAGACCG | ZNF470                  |              388566 |
        """
    },
    {
        "role": "assistant",
        "content": """
        {
          "article_title": "Cooperation of the ATM and Fanconi Anemia/BRCA Pathways in Double-Strand Break End Resection",
          "dataset_filename": "caiCooperationATMFanconi2020_32075772_NIHMS1563617-supplement-2_CP0041_20170705_compat_chip_t1.csv",
          "header_row": "0",
          "sgRNA_sequence": "None",
          "gene_identifier": "HGNC gene symbol",
          "dataset_type": "Lookup/Reference/Metadata",
          "statistic_aggregation": "NA",
          "metrics": ["NA"],
          "conditions": ["NA"]
        }
        """
    },
    {
        "role": "user",
        "content": """
        **Article:** A Whole-genome CRISPR Screen Identifies a Role of MSH2 in Cisplatin-mediated Cell Death in Muscle-invasive Bladder Cancer
        **sample:**
        | Supplementary Table 5. Raw counts from the cisplatin resistance screen.    | Unnamed: 1   | Unnamed: 2   | Unnamed: 3   | Unnamed: 4   | Unnamed: 5   | Unnamed: 6   | Unnamed: 7   | Unnamed: 8   |
        |:---------------------------------------------------------------------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|
        | sgRNA                                                                      | Cisplatin1   | Cisplatin2   | Cisplatin3   | Cisplatin4   | DMSO1        | DMSO2        | DMSO3        | DMSO4        |
        | EHMT1-HGLibA_14685                                                         | 123          | 228          | 147          | 16           | 224          | 234          | 208          | 89           |
        | RALGPS2-HGLibA_40172                                                       | 5            | 2            | 7            | 167          | 15           | 155          | 36           | 90           |
        | HIST2H2BE-HGLibA_21534                                                     | 9            | 492          | 1            | 49           | 28           | 0            | 72           | 91           |
        | SGCA-HGLibA_43795                                                          | 173          | 87           | 273          | 82           | 142          | 111          | 174          | 156          |
        """
    },
    {
        "role": "assistant",
        "content": """
        {
          "article_title": "A Whole-genome CRISPR Screen Identifies a Role of MSH2 in Cisplatin-mediated Cell Death in Muscle-invasive Bladder Cancer",
          "dataset_filename": "goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 5_t1.csv",
          "header_row": "1",
          "sgRNA_sequence": "None",
          "gene_identifier": "HGNC gene symbol [MOD]",
          "dataset_type": "Data",
          "statistic_aggregation": "Per sgRNA",
          "metrics": ["Raw count"],
          "conditions": ["Cisplatin", "DMSO"]
        }
        """
    },
]

def process_dataset(df, filename, paper_title, max_retries=5, output_file='gpt_output/output_4.md'):
    time.sleep(1)
    # Generate the first 10 lines of the dataset as string
    dataset_sample = df.head(10).to_string()

    # Generate API call
    user_input = [
          {"role": "system", "content": "You are an app that understands and analyzes CRISPR screen datasets and outputs structured JSON data specifying key parameters of the dataset in human readable form. All of the datasets you will see pertain to screens looking for genes which promote resistance (or sensitivity) to DNA damaging insults (drugs or radiation)."},
          {"role": "user", "content": f"I need to analyze a dataset from this paper: '{paper_title}'.\n\nThe filename of this dataset is '{filename}'. Here are the first 10 lines of the dataset:\n\n{dataset_sample}\n\nI need an interpretation of the dataset in structured JSON format. (A few examples of correcly labeled datasets are shown below.)\n\n"}
      ]
    
    message_input = user_input + prompt_examples
    
    retries = 0
    while retries < max_retries:
        try:
            response = openai.ChatCompletion.create(
                model="gpt-4-0613",
                messages = message_input,
                functions = [
                    {
                    "name": "analyze_dataset",
                    "description": "Analyzes a CRISPR screen dataset and provides structured data interpretation.",
                    "parameters": {
                      "type": "object",
                      "properties": {
                        "article_title": {
                          "type": "string",
                          "description": "The title of the paper"
                        },
                        "dataset_filename": {
                          "type": "string",
                          "description": "The filename of the dataset."
                        },
                        "header_row": {
                          "type": "string",
                          "description": "Which row the header is in (if any). The first row is '0'. If there is no header, then this should be 'None'."
                        },
                        "sgRNA_sequence": {
                          "type": "string",
                          "enum": [
                            "Own-column",
                            "Concatenated",
                            "None"
                          ],
                          "description": "Whether the sgRNA sequence is included in the table, and if so, whether on its own or concatenated with the gene name e.g. A1BG_CATCTTCTTTCACCTGAACG."
                        },
                        "gene_identifier": {
                          "type": "string",
                          "enum": [
                            "HGNC gene symbol",
                            "HGNC ID",
                            "Ensembl Gene ID",
                            "EntrezGene ID",
                            "RefSeq ID",
                            "UniProt ID",
                            "HGNC gene symbol [MOD]",
                            "HGNC ID [MOD]",
                            "Ensembl Gene ID [MOD]",
                            "EntrezGene ID [MOD]",
                            "RefSeq ID [MOD]",
                            "UniProt ID [MOD]",
                            "Other/Unclear [MOD]",
                            "None"
                          ],
                          "description": "The type of gene identifier used in the table. Add [MOD] if you recognize the gene identifier but it is modified in some way e.g. appended with a uuid (but not if it is simply appended with sgRNA sequence e.g. A1BG_CATCTTCTTTCACCTGAACG. In this case do not add [MOD])."
                        },
                        "dataset_type": {
                          "type": "string",
                          "enum": [
                            "Data",
                            "Lookup/Reference/Metadata"
                          ],
                          "description": "Whether the dataset contains experimentally measured data or is merely a lookup/reference/metadata table."
                        },
                        "statistic_aggregation": {
                          "type": "string",
                          "enum": [
                            "Per gene",
                            "Per sgRNA",
                            "NA"
                          ],
                          "description": "Whether the statistics are aggregated per gene or per sgRNA. If there is more than one row for any gene, then this should be 'per sgRNA', even if sgRNA sequence is not explicitly present. NA if the dataset is of type 'Lookup/Reference/Metadata'."
                        },
                        "metrics": {
                          "type": "array",
                          "items": {
                            "type": "string",
                            "enum": [
                              "MAGeCK statistics",
                              "Z-score/Z-ratios",
                              "fold-change/log fold-change",
                              "STARS (negative binomial)",
                              "Raw count",
                              "Other/ambiguous",
                              "NA"
                            ]
                          },
                          "description": "The metrics of gene/sgRNA abundance in the dataset (all that apply, but don't mention fold change if it's part of MAGeCK). NA if the dataset is of type 'Lookup/Reference/Metadata'."
                        },
                        "conditions": {
                          "type": "array",
                          "items": {
                            "type": "string"
                          },
                          "description": "The conditions explicitly mentioned in the dataset (or inferred from the file name if high confidence), e.g., DMSO, library, cisplatin, etc, if any. None if there are no explicitly mentioned conditions. NA if the dataset is of type 'Lookup/Reference/Metadata'."
                        }
                      }
                    },
                    "required": [
                      "article_title",
                      "dataset_filename",
                      "header_row",
                      "sgRNA_sequence",
                      "gene_identifier",
                      "dataset_type",
                      "statistic_aggregation",
                      "metrics",
                      "conditions"
                    ]
                  }
                ],
                function_call={"name": "analyze_dataset"},
                temperature=0
            )

            # Convert the response to a dictionary
            response_dict = json.loads(str(response.choices[0]["message"]["function_call"]["arguments"]))

            df.columns = df.columns.str.replace('|', '_')
            df.replace('\|', '_', regex=True, inplace=True, )
            dataset_sample = df.head(10).to_markdown()
            
            # Add the dataset_sample to the dictionary
            response_dict["sample"] = dataset_sample

            # Convert the updated dictionary 
            markdown_output = convert_to_markdown(response_dict)
            with open(output_file, 'a') as file:
                file.write(markdown_output)
                file.write("\n")
            break
        except Exception as e:
            print(f"Error: {e}")
            retries += 1
            if retries < max_retries:
                time.sleep(3)  # wait for 3 seconds before retrying
            else:
                print(
                    f"Failed to process dataset after {max_retries} attempts. Exiting.")
                break


# Load the BibTeX data and create a dictionary mapping pmid to title
with open('chemo_radio_resistance_screens.bib') as bibfile:
    bib_data = bibtexparser.load(bibfile)
    bib_dict = {entry['pmid']: entry['title'] for entry in bib_data.entries}

# Directory containing the datasets
dataset_directory = 'ifitm_datasets'

# Iterate over all CSV files in the directory
for filename in os.listdir(dataset_directory):
    if filename.endswith(".csv"):
        dataset_filepath = os.path.join(dataset_directory, filename)
        print(f"Processing dataset: {dataset_filepath}")

        # Extract pmid from filename
        pmid = filename.split('_')[1]

        # Find corresponding paper title from bib_dict
        if pmid in bib_dict:
            paper_title = bib_dict[pmid]

            # Analyze the dataset
            analyze_dataset(dataset_filepath, paper_title)
        else:
            print(f"No matching paper title found for pmid: {pmid}")
