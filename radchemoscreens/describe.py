import time
import os
import pandas as pd
import openai
import bibtexparser
import json

openai.api_key = "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"


def analyze_dataset(dataset_filepath, other_dataset_filenames, paper_title, paper_abstract):
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

    if file_ext == '.xlsx':
        for sheet_name, df in dfs.items():
            process_dataset(df, dataset_filepath, other_dataset_filenames,
                            paper_title, paper_abstract, sheet_name)
    else:
        process_dataset(df, dataset_filepath,
                        other_dataset_filenames, paper_title, paper_abstract, )


def process_dataset(df, dataset_filepath, other_dataset_filenames, paper_title, paper_abstract, sheet_name=None, max_retries=5, output_file='output.txt'):
    # Generate the first 20 lines of the dataset as string
    dataset_sample = df.head(30).to_string()

    # Generate API call
    
    if sheet_name:    
      message_input = [
          {"role": "system", "content": "You are an app that understands and analyzes CRISPR screen datasets and outputs structured JSON data specifying key parameters of the dataset in human readable form. All of the datasets you will see pertain to screens looking for genes which promote resistance (or sensitivity) to DNA damaging insults (drugs* or radiation).\n\n*drugs may include cisplatin, carboplatin, oxaliplatin, doxorubicin, daunorubicin, epirubicin, idarubicin, bleomycin, etoposide, irinotecan, topotecan, gemcitabine, 5-fluorouracil, capecitabine, methotrexate, pemetrexed, cyclophosphamide, ifosfamide, mitomycin C, camptothecin, methyl methanesulfonate, temozolomide, chlorambucil, busulfan, or other DNA damaging agents."},
          {"role": "user", "content": f"I need to analyze a dataset from this paper: '{paper_title}'. Here is its abstract: \n\n'{paper_abstract}'\n\nThe filename of this dataset is '{dataset_filepath}' - this is the sheet named '{sheet_name}', and there are also the following related datasets associated with this study: {other_dataset_filenames}. Here are the first 20 lines of the dataset:\n\n{dataset_sample}\n\nI need an interpretation of the dataset in structured JSON format. You can use the filename, abstract, and column names, as well as your background understanding of CRISPR screen methods, DNA damage, and data analysis, to make informed decisions about the relevant information. Be as informative as possible, but if there is anything that is particularly ambiguous, you can add an [unclear] tag to the end of the line."}
      ]
    else:
      message_input = [
        {"role": "system", "content": "You are an app that understands and analyzes CRISPR screen datasets and outputs structured JSON data specifying key parameters of the dataset in human readable form. All of the datasets you will see pertain to screens looking for genes which promote resistance (or sensitivity) to DNA damaging insults (drugs* or radiation).\n\n*drugs may include cisplatin, carboplatin, oxaliplatin, doxorubicin, daunorubicin, epirubicin, idarubicin, bleomycin, etoposide, irinotecan, topotecan, gemcitabine, 5-fluorouracil, capecitabine, methotrexate, pemetrexed, cyclophosphamide, ifosfamide, mitomycin C, camptothecin, methyl methanesulfonate, temozolomide, chlorambucil, busulfan, or other DNA damaging agents."},
        {"role": "user", "content": f"I need to analyze a dataset from this paper: '{paper_title}'. Here is its abstract: \n\n'{paper_abstract}'\n\nThe filename of this dataset is '{dataset_filepath}', and there are also the following related datasets associated with this study: {other_dataset_filenames}. Here are the first 20 lines of the dataset:\n\n{dataset_sample}\n\nI need an interpretation of the dataset in structured JSON format. You can use the filename, abstract, and column names, as well as your background understanding of CRISPR screen methods, DNA damage, and data analysis, to make informed decisions about the relevant information. Be as informative as possible, but if there is anything that is particularly ambiguous, you can add an [unclear] tag to the end of the line."}
      ]
        
    
    retries = 0
    while retries < max_retries:
        try:
            response = openai.ChatCompletion.create(
                model="gpt-3.5-turbo-16k-0613",
                messages=message_input,
                functions= [
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
                                "sheet_name": {
                                  "type": "string",
                                  "enum": [f"{sheet_name}"]
                                },
                                "dataset_type": {
                                    "type": "string",
                                    "enum": ["All genes/guides", "Hits only"],
                                    "description": "The type of dataset (pick one)"
                                },
                                "metrics": {
                                    "type": "array",
                                    "items": {
                                        "type": "string",
                                        "enum": ["P-value/FDR", "score", "fold change", "count"]
                                    },
                                    "description": "The metrics of gene/sgRNA abundance in the dataset (all that apply)"
                                },
                                "table_multiplicity": {
                                    "type": "string",
                                    "enum": ["Singular", "Multiple"],
                                    "description": "Whether multiple tables are being shown on a single sheet or not (e.g., if there are multiple 'empty' columns, and repeating column names, then there are multiple tables on a single sheet)"
                                },
                                "table_name_or_names": {
                                    "type": "array",
                                    "items": {"type": "string"},
                                    "description": "Proposed name (or list of names) for the table(s) in the sheet."
                                },
                                "col_names_interpretations": {
                                    "type": "array",
                                    "items": {
                                        "type": "object",
                                        "properties": {
                                            "dataset_filename_table_name_column_name": {
                                                "type": "string",
                                                "description": "The name of the dataset, appended to the name of the table, appended to the actual column name."
                                            },
                                            "interpretation": {
                                                "type": "string",
                                                "description": "Proposed interpretation of what the column represents."
                                            }
                                        }
                                    },
                                    "description": "A list of tables/dataframes with two columns: `dataset_filename_table_name_column_name` and `interpretation`."
                                }
                            },
                            "required": [
                                "article_title",
                                "dataset_filename",
                                "sheet_name",
                                "dataset_type",
                                "metrics",
                                "table_multiplicity",
                                "table_name_or_names",
                                "col_names_interpretations"
                            ]
                        }
                    }
                ],
                function_call={"name": "analyze_dataset"},
                temperature=0
            )

            with open(output_file, 'a') as file:
                file.write(str(response.choices[0]["message"]["function_call"]["arguments"]))
                file.write("\n")
            break
        except Exception as e:
            print(f"Error: {e}")
            retries += 1
            if retries < max_retries:
                time.sleep(2)  # wait for 2 seconds before retrying
            else:
                print(
                    f"Failed to process dataset after {max_retries} attempts. Exiting.")
                break


# Load the BibTeX data and iterate over the entries
with open('chemo_radio_resistance_screens.bib') as bibfile:
    bib_data = bibtexparser.load(bibfile)

for entry in bib_data.entries:
    bib_key = entry['ID']
    paper_title = entry['title']
    paper_abstract = entry['abstract']
    pmid = entry['pmid']

    # Search for related datasets
    dataset_directory = 'downloads'
    dataset_files = [os.path.join(dataset_directory, filename) for filename in os.listdir(dataset_directory) if (pmid in filename or bib_key in filename)]
    if dataset_files:
            for i, primary_dataset in enumerate(dataset_files):
                other_datasets = dataset_files[:i] + dataset_files[i+1:]
                analyze_dataset(primary_dataset, other_datasets, paper_title, paper_abstract)