import os
import sys
import argparse
import pandas as pd
from tabulate import tabulate

def main():

    # Define the argument parser
    parser = argparse.ArgumentParser(description='Generate previews for CSV files.')
    parser.add_argument('-d', '--dir', help='Directory containing the CSV files.')
    parser.add_argument('-b', '--base', help='Base directory for files from stdin or queryfile.')
    parser.add_argument('-f', '--files', nargs='*', help='List of CSV files.')
    parser.add_argument('-qf', '--queryfile', help='Query file with a list of CSV files.')
    args = parser.parse_args()

    base_directory = args.base if args.base else ''
    
    # Use either the directory, the list of files or the queryfile
    if args.dir:
        directory_path = args.dir
        files = [f for f in os.listdir(directory_path) if f.endswith(".csv")]
        files.sort()
    elif args.files:
        files = [os.path.join(base_directory, f) for f in args.files]
    elif args.queryfile:
        if args.queryfile == '-':
            files = [os.path.join(base_directory, line.strip()) for line in sys.stdin]
        else:
            with open(args.queryfile, 'r') as f:
                files = [os.path.join(base_directory, line.strip()) for line in f]
    else:
        print("Please provide a directory (-d), a list of files (-f), or a queryfile (-qf).")
        sys.exit(1)

    md = open("key_stats.md", "w")

    # Iterate through all CSVs
    for file in files:
        df = pd.read_csv(file)

        n_rows = df.shape[0]
        n_cols = df.shape[1]
        preview = tabulate(df.head(), headers='keys', tablefmt='pipe', showindex=False)

        md.write(f'## {file}\n')
        md.write(f'**n_rows:** {n_rows}\n')
        md.write(f'**n_cols:** {n_cols}\n')
        md.write('**preview:**\n')
        md.write(preview)
        md.write('\n\n')

    md.close()

if __name__ == "__main__":
    main()
