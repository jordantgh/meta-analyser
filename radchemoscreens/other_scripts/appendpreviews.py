import os
import json
import pandas as pd

# Your output file
output_file = 'output_2.json'

# Load the entries from the output file
with open(output_file, 'r') as f:
    entries = json.load(f)

# Processed entries
processed_entries = []

for entry_dict in entries:
    # Load the dataset
    dataset_filename = os.path.join('output_tables', entry_dict['dataset_filename'] + '.csv')
    df = pd.read_csv(dataset_filename)

    # Get the first 5 lines of the dataset
    dataset_sample = df.head(5)

    # Format the dataset sample in markdown table format
    markdown_table = dataset_sample.to_markdown()

    # Add the markdown_table to the entry dict
    entry_dict['preview'] = markdown_table

    # Add the processed entry to the list
    processed_entries.append(entry_dict)

# Your new output file
new_output_file = 'output_3.json'

# Write the processed entries to the new output file
with open(new_output_file, 'w') as f:
    json.dump(processed_entries, f, indent=4)
