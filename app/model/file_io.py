import pandas as pd
import requests
import os
import xlwings as xw
import numpy as np


def extract_dfs(fname: 'str', should_stop: 'bool') -> 'dict[str, pd.DataFrame]':
    data = {}
    if should_stop:
        return

    if fname.endswith(('.xlsx')):
        if should_stop:
            return
        with xw.Book(fname, mode="r") as book:
            for sheet in book.sheets:
                data[sheet.name] = sheet.cells.options(
                    pd.DataFrame, header=False, index=False
                ).value

    elif fname.endswith('.csv'):
        if should_stop:
            return
        data = {"sheet": pd.read_csv(fname, header=None, index_col=None)}
    else:
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
                if should_stop:
                    os.remove(local_fname)
                    return
                f.write(chunk)

    return local_fname
