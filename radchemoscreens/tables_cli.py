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

def find_blocks(mask):
    # Copy the matrix to not modify the original
    matrix = np.copy(mask)

    # Define the rectangles list
    rectangles = []

    # Function to create a rectangle
    def create_rectangle(i, j):
        # Determine the height of the rectangle
        height = 0
        while i + height < len(matrix) and matrix[i + height][j] == 1:
            height += 1

        # Determine the width of the rectangle
        width = 0
        while j + width < len(matrix[0]) and (
            all(matrix[k][j + width] == 1 for k in range(i, i + height)) and
            (i + height == len(matrix) or matrix[i + height][j + width] == 0) and
            (i == 0 or matrix[i - 1][j + width] == 0)
        ):
            width += 1

        # Remove the rectangle cells from the matrix
        for x in range(i, i + height):
            for y in range(j, j + width):
                matrix[x][y] = 0

        # Return the rectangle if it's at least 1x2 and connected to an edge
        if height >= 1 and width >= 2 and (
            i == 0 or j == 0 or i + height - 1 == len(matrix) - 1 or j + width - 1 == len(matrix[0]) - 1
        ):
            return ((i, j), (i + height - 1, j + width - 1))
        else:
            return None

    # Loop over every cell in the matrix
    for j in range(len(matrix[0])):
        for i in range(len(matrix)):
            # If the cell is a '1', create a rectangle
            if matrix[i][j] == 1:
                rectangle = create_rectangle(i, j)
                if rectangle is not None:
                    rectangles.append(rectangle)

    # Return the list of rectangles
    return rectangles

def parse_tables(df):
    binary_rep = np.array(df.notnull().astype("int"))
    list_of_dataframes = []
    labeled = label(binary_rep)

    for region in regionprops(labeled):
        minr, minc, maxr, maxc = region.bbox
        r = region # debug
        region = df.iloc[minr:maxr, minc:maxc]
        
        print("Region:", r.label)  # print the label of the region
        print(region.head(2))  # print the first two rows of the region dataframe


        # Binarise the region and invert (so empty cells are `1` and cells with values are `0`)
        binary_region = np.array(region.isnull().astype("int"))

        # Apply the find_blocks function to get the list of blocks
        blocks = find_blocks(binary_region)
        
        top_blocks = []
        bottom_blocks = []
        top_block_columns = []

        # Assign blocks to top_blocks and bottom_blocks and track column ranges for top blocks
        for block in blocks:
            if block[0][0] == 0:  # Top block starts from row 0
                top_blocks.append(block)
                top_block_columns.append(set(range(block[0][1], block[1][1] + 1)))
            elif block[1][0] == binary_region.shape[0] - 1:  # Bottom block ends at the last row
                bottom_blocks.append(block)
                
        if not top_blocks:
            list_of_dataframes.append(region)
            continue

        for block in top_blocks:
            # Create initial table using top block.
            top_rows_removed = block[1][0] - block[0][0] + 1  # total rows in top block
            rowrange = slice(block[1][0] + 1, binary_region.shape[0])
            colrange = slice(block[0][1], block[1][1] + 1)
            subregion = region.iloc[rowrange, colrange]

            # Check if a bottom block's column indices match exactly with the subregion's column indices.
            for bblock in bottom_blocks:
                bcolrange = slice(bblock[0][1], bblock[1][1] + 1)
                if colrange == bcolrange:
                    # If there is an exact match, slice the subregion further to exclude the bottom block.
                    # Adjust the bottom block slice to be relative to the top block slice.
                    browrange = bblock[0][0] - top_rows_removed
                    subregion = subregion.iloc[:browrange, :]

            list_of_dataframes.append(subregion)

        # After processing top blocks, find any column ranges not covered by top blocks and treat these as tables
        all_columns = set(range(binary_region.shape[1]))
        covered_columns = set.union(*top_block_columns)
        uncovered_columns = all_columns - covered_columns

        while uncovered_columns:
            start_col = min(uncovered_columns)
            end_col = start_col
            while end_col + 1 in uncovered_columns:
                end_col += 1

            uncovered_columns -= set(range(start_col, end_col + 1))
            colrange = slice(start_col, end_col + 1)
            subregion = region.iloc[:, colrange]

            # Check if a bottom block's column indices match exactly with the subregion's column indices.
            for bblock in bottom_blocks:
                bcolrange = slice(bblock[0][1], bblock[1][1] + 1)
                if colrange == bcolrange:
                    # If there is an exact match, slice the subregion further to exclude the bottom block.
                    # Adjust the bottom block slice to be relative to the top block slice.
                    browrange = bblock[0][0]
                    subregion = subregion.iloc[:browrange, :]

            list_of_dataframes.append(subregion)

    return list_of_dataframes

def main():
    parser = argparse.ArgumentParser(description="Parse and save tables from Excel files.")
    parser.add_argument("file_paths", nargs="+", help="One or more file paths for processing.")
    args = parser.parse_args()

    # Loop over each specified file
    for file_path in args.file_paths:
        if not os.path.isfile(file_path):
            print(f"Error: File {file_path} does not exist.")
            continue

        try:
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
                  
                    filename = f"{sheet_name}_df_{i}.csv"  # Generate a unique filename for each DataFrame
                    df.to_csv(filename, index=False, header=False)
                    print(f"Saved DataFrame {i} from sheet {sheet_name} as {filename}")

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")

if __name__ == "__main__":
    main()