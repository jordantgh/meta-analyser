import pickle

from model.article_managers import (
    Bibliography, Article, ProcessedTable,
    ProcessedTableManager, stash_all_observers, restore_all_observers
)

from model.database import TableDBManager, PostPruningTableDBEntry
from model.threading import (
    SearchThread, FilePreviewThread, FileProcessingThread
)

from utils.constants import PageIdentity, Mode
import scripts.query_parser as qp


class Model:
    def __init__(self):
        self._state = Mode.BROWSING
        self.bibliography = Bibliography()
        self.search_thread = SearchThread()
        self.search_preview_thread = FilePreviewThread("")
        self.table_db_manager = TableDBManager()
        self.processed_table_manager = ProcessedTableManager()
        self.processing_thread = FileProcessingThread(self.table_db_manager)

        # Flags to track whether we've ever parsed/pruned for when we load
        self.n_parse_runs = 0
        self.n_prunes = 0

    @property
    def state(self):
        return self._state

    def set_state(self, state):
        self._state = state

        if state == Mode.PROCESSING:
            self.n_parse_runs += 1
        elif state == Mode.PRUNING:
            self.n_prunes += 1

    def create_article_data(self, article_json):
        article = Article(article_json)
        self.bibliography.add_article(article)

        return article

    def _update_processed_tables(self, article, ids_list):
        processed_tables = []
        for table_id, file_id in ids_list:
            table_data = self.table_db_manager.get_processed_table_data(
                table_id,
                PageIdentity.PARSED
            )

            if table_data is not None:
                num_columns = len(table_data.columns)
            else:
                num_columns = None

            processed_table = ProcessedTable(
                article, table_id, file_id, num_columns
            )

            self.processed_table_manager.add_processed_table(processed_table)
            processed_tables.append(processed_table)

        return processed_tables

    def update_article(self, article, ids_list):
        article.processed_tables = self._update_processed_tables(
            article, ids_list
        )

        return article

    def prune_tables_and_columns(self, context):
        for article in self.bibliography.get_selected_articles(context):

            # TODO unspaghettify this
            tables_to_prune = [table
                               for table
                               in article.processed_tables if table.checked]

            for table in tables_to_prune:
                pruned_df = None

                data = self.table_db_manager.get_processed_table_data(
                    table.id, context
                )

                if context == PageIdentity.PARSED:
                    cols = table.checked_columns
                elif context == PageIdentity.PRUNED:
                    cols = table.pruned_columns

                if data is not None and cols is not None:
                    pruned_df = data.iloc[:, cols]

                if pruned_df is not None:
                    existing_table = self.table_db_manager.get_table_object(
                        PostPruningTableDBEntry,
                        table.id
                    )

                    if existing_table:
                        self.table_db_manager.update_table(
                            PostPruningTableDBEntry,
                            table.id,
                            pruned_df
                        )

                    else:
                        self.table_db_manager.save_table(
                            PostPruningTableDBEntry,
                            table.id,
                            table.file_id,
                            pruned_df
                        )

            article.pruned_tables = tables_to_prune

            for table in tables_to_prune:
                new_data = self.table_db_manager.get_processed_table_data(
                    table.id, context
                )

                if new_data is not None:
                    table.pruned_columns = list(range(len(new_data.columns)))

    # TODO #30 filtering is slow and needs its own thread; gui hangs up too long
    def filter_tables(self, query, context):
        for article in self.bibliography.get_selected_articles(
            PageIdentity.PARSED
        ):
            for processed_table in article.processed_tables:
                table_data = self.table_db_manager.get_processed_table_data(
                    processed_table.id,
                    context
                )

                processed_table.set_checked(bool(qp.search(
                    query,
                    [(processed_table.id, table_data.to_string())])),
                    PageIdentity.PARSED
                )

    def reset_for_searching(self):
        self.bibliography.reset()

        self.n_parse_runs = 0
        self.n_prunes = 0

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
            'processed_db_url': self.table_db_manager.processed_db_url,
            'post_pruning_db_url': self.table_db_manager.post_pruning_db_url,
            'processed_table_manager': self.processed_table_manager,
            'n_parse_runs': self.n_parse_runs,
            'n_prunes': self.n_prunes
        }

        with open(filename, 'wb') as f:
            pickle.dump(save_object, f)

        visited_objects.clear()
        restore_all_observers(self.bibliography, global_stash, visited_objects)

    def load(self, filename):
        with open(filename, 'rb') as f:
            save_object = pickle.load(f)

        self.bibliography = save_object['bibliography']
        self.table_db_manager = TableDBManager(
            save_object['processed_db_url'],
            save_object['post_pruning_db_url'])
        self.processed_table_manager = save_object['processed_table_manager']
        self.search_thread = SearchThread()
        self.search_preview_thread = FilePreviewThread("")
        self.processing_thread = FileProcessingThread(self.table_db_manager)

        self.n_parse_runs = save_object.get('n_parse_runs', False)
        self.n_prunes = save_object.get('n_prunes', False)
