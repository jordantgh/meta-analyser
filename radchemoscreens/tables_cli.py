from skimage.measure import label, regionprops
import pandas as pd
import numpy as np
import argparse
import os

def load_data(file_path):
    # Check file extension to determine appropriate parsing method
    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)
    elif file_path.endswith(".txt"):
        df = pd.read_csv(file_path, delimiter="\t")  # Adjust delimiter as needed
    elif file_path.endswith(".xlsx"):
        df = pd.read_excel(file_path, sheet_name=None, header=None)  # Returns a dictionary
    else:
        raise ValueError("File format not supported.")
    return df

def parse_tables(df):
    binary_rep = np.array(df.notnull().astype("int"))
    labeled = label(binary_rep)

    list_of_dataframes = []
    for region in regionprops(labeled):
        minr, minc, maxr, maxc = region.bbox
        if maxr - minr <= 1 or maxc - minc <= 1:
            continue
        region = df.iloc[minr:maxr, minc:maxc]
        list_of_dataframes.append(region)

    return list_of_dataframes

def main():
    parser = argparse.ArgumentParser(description="Parse and save tables from Excel files.")
    parser.add_argument("file_paths", nargs="+", help="One or more file paths for processing.")
    args = parser.parse_args()

    # Output directory path
    output_dir = "output_tables"

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Loop over each specified file
    for file_path in args.file_paths:
        if not os.path.isfile(file_path):
            print(f"Error: File {file_path} does not exist.")
            continue

        try:
            # Get parent file name without extension
            parent_filename = os.path.splitext(os.path.basename(file_path))[0]

            data = load_data(file_path)

            # Check if the loaded data is a DataFrame or a dictionary of DataFrames (multi-sheet Excel file)
            if isinstance(data, pd.DataFrame):
                data = {"sheet": data}

            all_tables = {}

            # Iterating over all sheets
            for sheet_name, df in data.items():
                tables = parse_tables(df)
                all_tables[sheet_name] = tables  # Save the list of tables for each sheet

            # Save each DataFrame to a separate CSV file
            for sheet_name, tables in all_tables.items():
                for i, df in enumerate(tables):
                  
                    filename = f"{parent_filename}_{sheet_name}_t{i+1}.csv"  # Generate a unique filename for each DataFrame
                    output_path = os.path.join(output_dir, filename)  # Full output path
                    df.to_csv(output_path, index=False, header=False)
                    print(f"Saved DataFrame {i} from sheet {sheet_name} as {output_path}")

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")

if __name__ == "__main__":
    main()