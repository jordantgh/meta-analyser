import os
import logging
import requests
import pandas as pd
import numpy as np
from uuid import uuid4
from PyQt5.QtCore import QThread, pyqtSignal

import scripts.query_parser as qp

from skimage.measure import label, regionprops
from search_for_papers import query_pmc

from sqlalchemy import create_engine, Column, Integer, String, BLOB
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

import pickle

def save_processed_df_to_db(db_manager, unique_id, df):
    serialized_df = pickle.dumps(df)
    db_manager.save_table(unique_id, serialized_df)

Base = declarative_base()

class ProcessedTable(Base):
    __tablename__ = 'processed_tables'
    
    id = Column(Integer, primary_key=True)
    file_id = Column(String)
    table_data = Column(BLOB)


class TableDatabaseManager:
    def __init__(self, db_url="sqlite:///tables.db"):
        self.engine = create_engine(db_url)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

    def save_table(self, file_id, table_data):
        with self.Session() as session:
            new_table = ProcessedTable(file_id=file_id, table_data=table_data)
            session.add(new_table)
            session.commit()
            return new_table.id

    def get_table(self, file_id):
        with self.Session() as session:
            table = session.query(ProcessedTable).filter_by(file_id=file_id).first()
            if table:
                table = pickle.loads(table.table_data)
                table.reset_index(drop=True, inplace=True)
                return table
            return None

    def delete_table(self, table_id):
        with self.Session() as session:
            table = session.query(ProcessedTable).filter_by(id=table_id).first()
            if table:
                session.delete(table)
                session.commit()


def parse_tables(selected_articles, db_manager, callback=None):
    def _is_contained(bbox1, bbox2):
        return bbox2[0] <= bbox1[0] and bbox2[1] <= bbox1[1] and bbox2[2] >= bbox1[2] and bbox2[3] >= bbox1[3]

    for index,  article in enumerate(selected_articles):
        processed_table_ids = []
        for file in article.supp_files:
            if not file.checked: return
            
            try:
                fname = download_file(file.url)
                if fname is None: continue
                
                data = extract_data_from_file(fname)
                    
                for sheetname, df in data.items():
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
                        unique_id = f"{os.path.splitext(fname)[0]}_{sheetname}_Table{i}"
                        save_processed_df_to_db(db_manager, unique_id, region)
                        processed_table_ids.append(unique_id)

            except Exception as e:
                logging.exception(f"An error occurred while processing {file.url}: {str(e)}")

        if callback:
            progress = int(100 * (index + 1) / len(selected_articles))
            callback(article, processed_table_ids, progress)


def extract_data_from_file(fname): # returns a dictionary of dataframes
    if fname.endswith(('.xlsx')):
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
        query_pmc(self.query, callback=self.article_sig.emit, thread = self)
        self.finished_sig.emit()


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
    article_sig = pyqtSignal(object, list, int)
    finished_sig = pyqtSignal()

    def __init__(self, db_manager):
        super().__init__()
        self.selected_articles = []
        self.db_manager = db_manager

    def run(self):
        parse_tables(self.selected_articles, self.db_manager, callback=self.article_sig.emit)
        self.finished_sig.emit()


class SupplementaryFile:
    def __init__(self, article_id, url, id):
        self.checked = True
        self.article_id = article_id
        self.url = url
        self.id = id


class ProcessedTableFile:
    def __init__(self, article_id, id):
        self.checked = True
        self.article_id = article_id
        self.id = id
        self.checked_columns = set()
        
    def toggle_column_check(self, col_index):
        if col_index in self.checked_columns:
            self.checked_columns.remove(col_index)
        else:
            self.checked_columns.add(col_index)

class SupplementaryFileManager:
    def __init__(self):
        self.supp_files = {}

    def add_file(self, file):
        self.supp_files[file.id] = file

    def get_file(self, file_id):
        return self.supp_files.get(file_id)


class Article:
    def __init__(self, title, abstract, pmc_id, supp_files=[], tables=[]):
        self.checked = True
        self.title = title
        self.abstract = abstract
        self.pmc_id = pmc_id
        self.supp_files = supp_files
        self.processed_tables = tables

    def get_file(self, file_id):
        return next((f for f in self.supp_files if f.id == file_id), None)


class Bibliography:
    def __init__(self):
        self.articles = {}

    def add_article(self, article):
        self.articles[article.pmc_id] = article

    def get_article(self, article_id):
        return self.articles.get(article_id)

    def get_selected_articles(self):
        return [p for p in self.articles.values() if p.checked]


class Model:
    def __init__(self):
        self.filtered_articles = {}
        self.processing_mode = False
        self.bibliography = Bibliography()
        self.file_manager = SupplementaryFileManager()
        self.search_thread = SearchThread()
        self.preview_thread = FilePreviewThread("")   
        self.table_db_manager = TableDatabaseManager()
        self.processing_thread = FileProcessingThread(self.table_db_manager)

    def update_supp_files(self, article, article_json):
        supp_files = []
        for file_url in article_json["SupplementaryFiles"]:
            supp_file = SupplementaryFile(article.pmc_id, file_url, uuid4())
            self.file_manager.add_file(supp_file)
            supp_files.append(supp_file)
        article.supp_files = supp_files

    def create_article_data(self, article_json):
        article = Article(
            title=article_json["Title"],
            abstract=article_json["Abstract"],
            pmc_id=article_json["PMCID"]
        )

        self.update_supp_files(article, article_json)
        self.bibliography.add_article(article)

        return article

    def update_processed_tables(self, article, ids_list):
        processed_tables = []
        for table_id in ids_list:
            processed_table = ProcessedTableFile(article.pmc_id, table_id)
            processed_tables.append(processed_table)

        return processed_tables

    def update_article(self, article, ids_list):
        article.processed_tables = self.update_processed_tables(article, ids_list)

        return article

    def add_processed_tables(self, file_id, tables):
        file = self.file_manager.get_file(file_id)
        file.processed_tables = tables

    def filter_tables(self, query):
        self.filtered_articles = {}
        for article in self.bibliography.get_selected_articles():
            filtered_tables = []
            for processed_table in article.processed_tables:
                table_data = self.table_db_manager.get_table(processed_table.id)
                if qp.search(query, [(processed_table.id, table_data.to_string())]):
                    filtered_tables.append(processed_table)
            
            if filtered_tables:
                self.filtered_articles[article.pmc_id] = filtered_tables

    def apply_filtered_articles(self):
        for article_id, filtered_tables in self.filtered_articles.items():
            article = self.bibliography.get_article(article_id)
            if article:
                article.processed_tables = filtered_tables

