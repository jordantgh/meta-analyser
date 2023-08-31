#!/bin/bash

# Get script directory
script_dir=$(dirname "$0")

# Set your downloads directory
downloads_dir="$script_dir/downloads"

# Loop over every file in the downloads directory
for file in "$downloads_dir"/*
do
    # Execute the python script with the current file as argument
    python "$script_dir/tables_cli.py" "$file"
done