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

with open('dependencies/prompt_examples.json', 'r') as file:
    prompt_examples = json.load(file)
    
with open('dependencies/functioncall.json', 'r') as file:
    functioncall = json.load(file)

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
                functions = functioncall,
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
