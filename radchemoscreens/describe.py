import time
import os
import pandas as pd
import openai
import bibtexparser
import json
from dotenv import load_dotenv

load_dotenv()
openai.api_key = os.environ.get("OPENAI_API_KEY")

def analyze_dataset(dataset_filepath, paper_title):
    # Load the dataset
    file_ext = os.path.splitext(dataset_filepath)[-1]
    if file_ext == '.csv':
        df = pd.read_csv(dataset_filepath, header=None)
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
  
sys_prompt = "You're an app that analyzes CRISPR screen datasets and outputs structured JSON data with key parameters of the dataset. All datasets you'll see pertain to screens for genes promoting resistance or sensitivity to DNA damaging insults (drugs/radiation)."

with open('dependencies/prompt_examples.json', 'r') as file:
    cot = json.load(file)
    
with open('dependencies/functioncall.json', 'r') as file:
    functioncall = json.load(file)

def process_dataset(df, filename, paper_title, max_retries=5, output_file='gpt_output/output_5.md'):
    time.sleep(1)

    # Generate API call
    dataset_sample = df.head(6).to_markdown()
    user_prompt = [f"""
    I need to analyze a dataset from this paper: '{paper_title}'. The dataset filename is '{filename}'. Here are the first 5 rows:
    
    {dataset_sample}
    
    I need an interpretation of this dataset in structured JSON format. (A few examples of correcly labeled datasets are shown above.)
    """]

    prompt = [{"role": "system", "content": sys_prompt}, {"role": "user", "content": str(cot + user_prompt)}]
    
    message_input = prompt 
    
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

            # df.columns = df.columns.str.replace('|', '_')
            df.replace('\|', '_', regex=True, inplace=True, )
            dataset_sample = df.head(6).to_markdown()
            
            # Add the dataset_sample to the dictionary
            response_dict["sample"] = dataset_sample

            # Convert the updated dictionary 
            markdown_output = convert_to_markdown(response_dict)
            with open(output_file, 'a') as file:
                file.write(markdown_output)
                file.write("\n")
                
            # Write the prompt to a file
            # with open('gpt_output/prompt.md', 'w') as file:
            #     file.write(str(message_input))
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
