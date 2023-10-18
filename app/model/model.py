from uuid import uuid4
import pickle
from enum import Enum, auto

from model.article_managers import Bibliography, Article, SuppFile, SuppFileManager, ProcessedTable, ProcessedTableManager, stash_all_observers, restore_all_observers
from model.database import TableDBManager, PostPruningTableDBEntry
from model.threading import SearchThread, FilePreviewThread, FileProcessingThread
import scripts.query_parser as qp


class Model:
    class Mode(Enum):
        BROWSING = 0
        SEARCHING = auto()
        PROCESSING = auto()
        PRUNING = auto()

    def __init__(self):
        self._state = Model.Mode.BROWSING
        self.bibliography = Bibliography()
        self.file_manager = SuppFileManager()
        self.search_thread = SearchThread()
        self.preview_thread = FilePreviewThread("")
        self.table_db_manager = TableDBManager()
        self.processed_table_manager = ProcessedTableManager()
        self.processing_thread = FileProcessingThread(self.table_db_manager)
        
        # Flags to keep track of whether we've ever parsed/pruned for when we
        # load a model from a file; need to know to decide whether to populate
        # the parsed/pruned pages
        self.ever_parsed = 0
        self.ever_pruned = 0


    @property
    def state(self):
        return self._state
    
    def set_state(self, state):
        self._state = state
        
        # TODO: Reset this when we start a new search
        if state == Model.Mode.PROCESSING:
            self.ever_parsed += 1
        elif state == Model.Mode.PRUNING:
            self.ever_pruned += 1
    
    def update_supp_files(self, article, article_json):
        supp_files = []
        for file_url in article_json["SupplementaryFiles"]:
            supp_file = SuppFile(article, file_url, uuid4())
            self.file_manager.add_file(supp_file)
            supp_files.append(supp_file)

        article.supp_files = supp_files

    def create_article_data(self, article_json):
        article = Article(
            title=article_json["Title"],
            authors=article_json["Authors"],
            abstract=article_json["Abstract"],
            pmc_id=article_json["PMCID"],
            url=article_json["URL"])

        self.update_supp_files(article, article_json)
        self.bibliography.add_article(article)

        return article

    def update_processed_tables(self, article, ids_list):
        processed_tables = []
        for table_id, file_id in ids_list:
            table_data = self.table_db_manager.get_processed_table_data(
                table_id)

            if table_data is not None:
                num_columns = len(table_data.columns)
            else:
                num_columns = None
            processed_table = ProcessedTable(
                article,
                table_id,
                file_id,
                num_columns)

            self.processed_table_manager.add_processed_table(processed_table)
            processed_tables.append(processed_table)
        return processed_tables

    def update_article(self, article, ids_list):
        article.processed_tables = self.update_processed_tables(
            article, ids_list)

        return article

    def add_processed_tables(self, file_id, tables):
        file = self.file_manager.get_file(file_id)
        file.processed_tables = tables

    def prune_tables_and_columns(self, context):
        for article in self.bibliography.get_selected_articles(context):

            # TODO unspaghettify this
            tables_to_prune = [table
                               for table
                               in article.processed_tables if table.checked]

            for table in tables_to_prune:
                pruned_df = None

                if context == 'parsed':
                    columns_vector = table.checked_columns
                    data = self.table_db_manager.get_processed_table_data(
                        table.id)
                elif context == 'pruned':
                    columns_vector = table.pruned_columns
                    data = self.table_db_manager.get_post_pruning_table_data(
                        table.id)

                if data is not None and columns_vector is not None:
                    pruned_df = data.iloc[:, columns_vector]

                if pruned_df is not None:
                    existing_table = self.table_db_manager.get_table_object(
                        PostPruningTableDBEntry,
                        table.id)

                    if existing_table:
                        self.table_db_manager.update_table(
                            PostPruningTableDBEntry,
                            table.id,
                            pruned_df)
                    else:
                        self.table_db_manager.save_table(
                            PostPruningTableDBEntry,
                            table.id,
                            table.file_id,
                            pruned_df)

            article.pruned_tables = tables_to_prune

            for table in tables_to_prune:
                latest_data = self.table_db_manager \
                    .get_post_pruning_table_data(table.id)

                if latest_data is not None:
                    table.pruned_columns = list(range(len(
                        latest_data.columns)))

    def filter_tables(self, query):
        for article in self.bibliography.get_selected_articles('parsed'):
            for processed_table in article.processed_tables:
                table_data = self.table_db_manager.get_processed_table_data(
                    processed_table.id)

                processed_table.set_checked(bool(qp.search(
                    query,
                    [(processed_table.id, table_data.to_string())])),
                    'parsed')

    def reset_for_searching(self):
        # TODO / BUG currently seems to be a problem where previously processed
        # articles are not being reset properly, or at least we don't see the
        # tables anymore when we click on them in the GUI after a second round
        # of searching/processing
        self.bibliography.reset()
        self.file_manager.reset()

    def reset_for_processing(self):
        self.processed_table_manager.reset()
        self.table_db_manager.reset()


    def save(self, filename):
        # stash observers so we can pickle the model
        global_stash = {}
        # since articles refer to suppfiles/tables and vice versa, we need to
        # keep track of visited objects to avoid infinite recursion
        visited_objects = set()
        
        stash_all_observers(self.bibliography, global_stash, visited_objects)
        save_object = {
            'state': self.state,
            'bibliography': self.bibliography,
            'file_manager': self.file_manager,
            'table_db_manager_processed_db_url': self.table_db_manager.processed_db_url,
            'table_db_manager_post_pruning_db_url': self.table_db_manager.post_pruning_db_url,
            'processed_table_manager': self.processed_table_manager,
            'ever_parsed': self.ever_parsed,
            'ever_pruned': self.ever_pruned
        }

        # Serialize the dictionary and save it to a file
        with open(filename, 'wb') as f:
            pickle.dump(save_object, f)

        visited_objects.clear()
        restore_all_observers(self.bibliography, global_stash, visited_objects)
            
    def load(self, filename):
        # Load File
        with open(filename, 'rb') as f:
            save_object = pickle.load(f)

        self.bibliography = save_object['bibliography']
        self.file_manager = save_object['file_manager']
        self.table_db_manager = TableDBManager(
            save_object['table_db_manager_processed_db_url'],
            save_object['table_db_manager_post_pruning_db_url'])
        self.processed_table_manager = save_object['processed_table_manager']
        self.search_thread = SearchThread()
        self.preview_thread = FilePreviewThread("")
        self.processing_thread = FileProcessingThread(self.table_db_manager)
        
        self.ever_parsed = save_object.get('ever_parsed', False)
        self.ever_pruned = save_object.get('ever_pruned', False)