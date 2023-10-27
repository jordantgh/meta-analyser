import pandas as pd
import requests
import os


def extract_dfs(fname: 'str', should_stop: 'bool') -> 'dict[str, pd.DataFrame]':
    if should_stop:
        return

    if fname.endswith(('.xlsx')):
        # Check should_stop before a potentially long operation
        if should_stop:
            return
        data = pd.read_excel(fname, sheet_name=None,
                             header=None, index_col=None)
    elif fname.endswith('.csv'):
        # Check should_stop before a potentially long operation
        if should_stop:
            return
        data = {"sheet": pd.read_csv(fname, header=None, index_col=None)}
    else:
        # Check should_stop before a potentially long operation
        if should_stop:
            return
        data = {"sheet": pd.read_csv(
            fname, delimiter='\t', header=None, index_col=None)}

    if should_stop:
        return

    os.remove(fname)

    return data


def download_supp(url: 'str', should_stop: 'str') -> 'str':
    if should_stop:
        return

    local_fname = url.split('/')[-1]
    if not local_fname.endswith(('.txt', '.csv', '.xlsx')):
        return

    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_fname, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                # Check should_stop within the loop
                if should_stop:
                    os.remove(local_fname)
                    return
                f.write(chunk)

    return local_fname
