import os
import logging
import requests
import pandas as pd
import numpy as np
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from skimage.measure import label, regionprops
from search_for_papers import query_pmc


def parse_tables(data):
    def _is_contained(bbox1, bbox2):
        return bbox2[0] <= bbox1[0] and bbox2[1] <= bbox1[1] and bbox2[2] >= bbox1[2] and bbox2[3] >= bbox1[3]

    list_of_dataframes = []
    
    for df in data.values():
        if df.empty: continue
        binary_rep = np.array(df.notnull().astype("int"))
        labeled = label(binary_rep)
        region_bboxes = [region.bbox for region in regionprops(labeled)]

        for i, bbox1 in enumerate(region_bboxes):
            minr1, minc1, maxr1, maxc1 = bbox1
            if maxr1 - minr1 <= 1 or maxc1 - minc1 <= 1:
                continue
            if any(_is_contained(bbox1, bbox2) for j, bbox2 in enumerate(region_bboxes) if i != j):
                continue
            region = df.iloc[minr1:maxr1, minc1:maxc1]
            list_of_dataframes.append(region)

    return list_of_dataframes


def extract_data_from_file(fname):
    if fname.endswith(('.xlsx', '.xls')):
        data = pd.read_excel(fname, sheet_name=None, header=None, index_col=None)
    elif fname.endswith('.csv'):
        data = {"sheet": pd.read_csv(fname, header=None, index_col=None)}
    else:
        data = {"sheet": pd.read_csv(fname, delimiter='\t', header=None, index_col=None)}
    
    os.remove(fname)

    return data


def download_file(url):
    local_fname = url.split('/')[-1]
    if not local_fname.endswith(('.txt', '.csv', '.xlsx')):
        return
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_fname, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    return local_fname


class SearchThread(QThread):
    article_sig = pyqtSignal(dict, int)
    finished_sig = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.query = ""
        self.should_stop = False

    def stop(self):
        self.should_stop = True

    def run(self):
        try:
            query_pmc(self.query, callback=self.article_sig.emit, thread = self)
            self.finished_sig.emit()
        except Exception as e:
            logging.error(f"Unhandled exception in SearchThread: {e}")


class FilePreviewThread(QThread):
    prev_ready_sig = pyqtSignal(dict)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url

    def run(self):
        fname = download_file(self.file_url)
        data = extract_data_from_file(fname)
        self.prev_ready_sig.emit(data)


class FileProcessingThread(QThread):
    tables_ready_sig = pyqtSignal(dict)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url

    def run(self):
        fname = download_file(self.file_url)
        data = extract_data_from_file(fname)
        tables = parse_tables(data)
        saved_tables = self.write_tables_to_files(tables)
        self.tables_ready_sig.emit(saved_tables)

    def write_tables_to_files(self, tables):
        # Save each table as a new CSV and return their file names
        saved_files = {}
        for idx, table in enumerate(tables):
            filename = f"table_{idx}.csv"
            table.to_csv(filename)
            saved_files[filename] = table

        return saved_files


class SupplementaryFile:
    def __init__(self, url):
        self.url = url
        self.checked = True


class Paper:
    def __init__(self, title, abstract, pmc_id, files=[]):
        self.title = title
        self.abstract = abstract
        self.pmc_id = pmc_id
        self.files = [SupplementaryFile(url) for url in files]
        self.checked = True


class Bibliography:
    def __init__(self):
        self.papers = []

    def add_paper(self, paper):
        self.papers.append(paper)
        
    def get_selected_papers(self):
        return [p for p in self.papers if p.checked]


class Model:
    def __init__(self):
        self.bibliography = Bibliography()
        self.search_thread = SearchThread()
        self.preview_thread = FilePreviewThread("")   
        
    def create_paper_data(self, article):
        return Paper(
            title=article["Title"],
            abstract=article["Abstract"],
            pmc_id=article["PMCID"],
            files=article["SupplementaryFiles"]
        )

    def extract_paper_details(self, item):
        paper = item.data(Qt.UserRole)
        paper_details = {
            'title': paper.title,
            'abstract': paper.abstract,
            'files': paper.files
        }

        return paper_details

    def extract_selected_papers_and_files(self):
        selected_papers = {}
        for paper in self.bibliography.get_selected_papers():
            selected_files = [f.url for f in paper.files if f.checked]
            if selected_files:
                selected_papers[paper.pmc_id] = {
                    "title": paper.title,
                    "files": selected_files
                }

        return selected_papers