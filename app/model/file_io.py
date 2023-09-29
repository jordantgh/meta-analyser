import pandas as pd
import requests
import os

def extract_dfs(fname): # returns a dictionary of dataframes
    if fname.endswith(('.xlsx')):
        data = pd.read_excel(fname, sheet_name=None, header=None, index_col=None)
    elif fname.endswith('.csv'):
        data = {"sheet": pd.read_csv(fname, header=None, index_col=None)}
    else:
        data = {"sheet": pd.read_csv(fname, delimiter='\t', header=None, index_col=None)}
    
    os.remove(fname)

    return data


def download_supp(url):
    local_fname = url.split('/')[-1]
    if not local_fname.endswith(('.txt', '.csv', '.xlsx')):
        return
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_fname, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    return local_fname
