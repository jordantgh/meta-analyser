from uuid import uuid4
import pickle

from model.article_managers import Bibliography, Article, SuppFile, SuppFileManager, ProcessedTable, ProcessedTableManager
from model.database import TableDBManager, PostPruningTableDBEntry
from model.threading import SearchThread, FilePreviewThread, FileProcessingThread
import scripts.query_parser as qp

class Model:
    def __init__(self):
        self.article_list_filtered = {}
        self.processing_mode = False
        self.bibliography = Bibliography()
        self.file_manager = SuppFileManager()
        self.search_thread = SearchThread()
        self.preview_thread = FilePreviewThread("")   
        self.table_db_manager = TableDBManager()
        self.processed_table_manager = ProcessedTableManager()
        self.processing_thread = FileProcessingThread(self.table_db_manager)

    def update_supp_files(self, article, article_json):
        supp_files = []
        for file_url in article_json["SupplementaryFiles"]:
            supp_file = SuppFile(article.pmc_id, file_url, uuid4())
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
        for table_id, file_id in ids_list:
            table_data = self.table_db_manager.get_processed_table_data(table_id)
            if table_data is not None:
                num_columns = len(table_data.columns)
            else:
                num_columns = None
            processed_table = ProcessedTable(article.pmc_id, table_id, file_id, num_columns)
            self.processed_table_manager.add_processed_table(processed_table)
            processed_tables.append(processed_table)
        return processed_tables

    def update_article(self, article, ids_list):
        article.processed_tables = self.update_processed_tables(article, ids_list)

        return article

    def add_processed_tables(self, file_id, tables):
        file = self.file_manager.get_file(file_id)
        file.processed_tables = tables

    def prune_tables_and_columns(self):
        for article in self.bibliography.get_selected_articles():
            ids = [table.id for table in self.article_list_filtered.get(article.pmc_id, []) if table.checked]
            
            tables_to_prune = []

            if article.to_prune:
                tables_to_prune = [self.processed_table_manager.get_processed_table(table_id) for table_id in article.to_prune if self.processed_table_manager.get_processed_table(table_id).checked]
            else:
                tables_to_prune = [table for table in article.processed_tables if table.checked]

            for table in tables_to_prune:
                serialized_df = None
                table_id = table.id
                data = self.table_db_manager.get_processed_table_data(table_id)

                if data is not None and table.checked_columns is not None:
                    data = data.iloc[:, table.checked_columns]
                    serialized_df = pickle.dumps(data)

                if serialized_df is not None:
                    existing_table = self.table_db_manager.get_table_object(PostPruningTableDBEntry, table_id)
                    if existing_table:
                        self.table_db_manager.update_table(PostPruningTableDBEntry, table_id, serialized_df)
                    else:
                        self.table_db_manager.save_table(PostPruningTableDBEntry, table_id, table.file_id, serialized_df)

        
            kept_tables = [self.processed_table_manager.get_processed_table(id) for id in ids if self.processed_table_manager.get_processed_table(id) is not None]

            for table in kept_tables:
                if table.checked_columns is not None:
                    latest_data = self.table_db_manager.get_processed_table_data(table.id)
                    if latest_data is not None:
                        table.checked_columns = list(range(len(latest_data.columns)))
              
    def filter_tables(self, query): 
        for article in self.bibliography.get_selected_articles():
            matched_tables = []
            for processed_table in article.processed_tables:
                table_data = self.table_db_manager.get_processed_table_data(processed_table.id)
                
                if table_data is not None:
                    if qp.search(query, [(processed_table.id, table_data.to_string())]):
                        matched_tables.append(processed_table)
                
            if matched_tables:
                self.article_list_filtered[article.pmc_id] = matched_tables

    def update_article_prune_list(self):
        for article_id, matched_tables in self.article_list_filtered.items():
            article = self.bibliography.get_article(article_id)
            if article:
                article.to_prune = matched_tables

    def reset_for_searching(self):
        self.bibliography.reset()
        self.file_manager.reset()

    def reset_for_processing(self):
        self.article_list_filtered = {}
        self.processed_table_manager.reset()
        self.table_db_manager.reset()