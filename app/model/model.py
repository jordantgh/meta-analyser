from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from uuid import UUID
    from pandas import DataFrame

import pickle
import os
from model.article_managers import (
    Bibliography, Article, ProcessedTable,
    ProcessedTableManager, stash_all_observers, restore_all_observers
)

from model.database import TableDBManager, PostPruningTableDBEntry, SortedListTableDBEntry
from model.threading import (
    SearchThread, FilePreviewThread, FileProcessingThread
)

from utils.constants import PageIdentity, Mode
import scripts.query_parser as qp
from model.tabular_operations import cols_to_sorted_lists


class Model:
    def __init__(
            self,
            db_temp_path_root: 'str' = None,
            db_perm_path_root: 'str' = None,
            saves_path: 'str' = None
    ):
        self.session_file = None
        self.db_temp_path_root = db_temp_path_root
        self.db_perm_path_root = db_perm_path_root
        self.db_session_files = None
        self.saves_path = saves_path
        self._state = Mode.BROWSING
        self.bibliography = Bibliography()
        self.search_thread = SearchThread()
        self.search_preview_thread = FilePreviewThread()
        self.table_db_manager = TableDBManager(
            db_temp_path_root, db_perm_path_root
        )
        self.processed_table_manager = ProcessedTableManager()
        self.processing_thread = FileProcessingThread(self.table_db_manager)
        self.last_selected_table: 'ProcessedTable' = None
        self.n_parse_runs = 0
        self.n_prunes = 0

    @property
    def state(self) -> 'Mode':
        return self._state

    def set_state(self, state: 'Mode'):
        self._state = state

        if state == Mode.PROCESSING:
            self.n_parse_runs += 1
        elif state == Mode.PRUNING:
            self.n_prunes += 1

    def create_article_data(self, article_json: 'dict'):
        article = Article(article_json)
        self.bibliography.add_article(article)

        return article

    def _update_processed_tables(
        self, article, ids_list: 'list[tuple[str, UUID]]'
    ):
        processed_tables: 'list[ProcessedTable]' = []
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

    def update_article(
        self, article: 'Article', ids_list: 'list[tuple[str, UUID]]'
    ):
        article.processed_tables = self._update_processed_tables(
            article, ids_list
        )

        return article

    # May need a thread for this, hangs up on large queries
    def prune_tables_and_columns(self, context: 'PageIdentity'):

        # TODO unspaghettify this
        article: 'Article'
        for article in self.bibliography.articles.values():

            unchecked_tables = [
                table for table in article.processed_tables if not table.checked
            ]

            for table in unchecked_tables:
                existing_table = self.table_db_manager.get_table_object(
                    PostPruningTableDBEntry, table.id
                )
                if existing_table:
                    self.table_db_manager.delete_table(
                        PostPruningTableDBEntry, table.id
                    )

            selected_tables = [
                table for table in article.processed_tables if table.checked
            ]

            table: 'ProcessedTable'
            for table in selected_tables:
                pruned_df = None

                data: 'DataFrame' = self.table_db_manager. \
                    get_processed_table_data(table.id, context)

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
                            article.pmc_id,
                            article.title,
                            article.url,
                            table.supp_file.url,
                            pruned_df,
                            table.id,
                            table.tags
                        )

            article.pruned_tables = selected_tables

            for table in article.pruned_tables:
                new_data = self.table_db_manager.get_processed_table_data(
                    table.id, PageIdentity.PRUNED
                )

                if new_data is not None:
                    table.pruned_columns = list(range(len(new_data.columns)))

    # TODO #30 filtering is slow and needs its own thread; gui hangs up too long
    def filter_tables(self, query: 'str', context: 'PageIdentity'):

        article: 'Article'
        for article in self.bibliography.get_selected_articles(
            PageIdentity.PARSED
        ):
            for processed_table in article.processed_tables:
                table_data = self.table_db_manager.get_processed_table_data(
                    processed_table.id,
                    context
                )

                processed_table.set_checked_state(
                    bool(
                        qp.search(
                            query, [
                                (processed_table.id, table_data.to_string())]
                        )),
                    PageIdentity.PARSED
                )

    def generate_sorted_lists(self):
        article: 'Article'
        for article in self.bibliography.get_selected_articles(
            PageIdentity.PRUNED
        ):
            for processed_table in article.processed_tables:
                if processed_table.mappings:
                    table_data = self.table_db_manager.get_processed_table_data(
                        processed_table.id,
                        PageIdentity.PRUNED
                    )

                    if table_data is not None:
                        sorted_lists = cols_to_sorted_lists(
                            processed_table, table_data
                        )
                    else:
                        continue

                    self.table_db_manager.save_table(
                        SortedListTableDBEntry,
                        article.pmc_id,
                        article.title,
                        article.url,
                        processed_table.supp_file.url,
                        sorted_lists,
                        processed_table.id,
                        processed_table.tags
                    )

    def reset_for_searching(self):
        self.bibliography.reset()

        self.n_parse_runs = 0
        self.n_prunes = 0

    def reset_for_processing(self):
        self.processed_table_manager.reset()
        self.table_db_manager.reset()

    def save(self, filename: 'str'):
        # stash observers so we can pickle the model
        global_stash = {}
        # since articles refer to suppfiles/tables and vice versa, we need to
        # keep track of visited objects to avoid infinite recursion
        visited_objects = set()

        stash_all_observers(self.bibliography, global_stash, visited_objects)

        if self.db_session_files:
            proc, prune, sorted_list = self.table_db_manager.save_dbs(
                filename, self.db_session_files
            )
        else:
            proc, prune, sorted_list = self.table_db_manager.save_dbs(filename)

        save_object = {
            'db_perm_path_root': self.db_perm_path_root,
            'state': self.state,
            'bibliography': self.bibliography,
            'processed_db_path': proc,
            'pruned_db_path': prune,
            'sorted_list_db_path': sorted_list,
            'processed_table_manager': self.processed_table_manager,
            'n_parse_runs': self.n_parse_runs,
            'n_prunes': self.n_prunes
        }

        with open(filename, 'wb') as f:
            pickle.dump(save_object, f)

        visited_objects.clear()
        restore_all_observers(self.bibliography, global_stash, visited_objects)
        self.session_file = filename
        self.db_session_files = [proc, prune, sorted_list]

    def load(self, file: 'str'):
        self.session_file = file
        self.saves_path = os.path.dirname(file)
        with open(file, 'rb') as f:
            save_object = pickle.load(f)
        self.bibliography = save_object['bibliography']
        self.db_session_files = [
            save_object['processed_db_path'],
            save_object['pruned_db_path'],
            save_object['sorted_list_db_path']
        ]
        self.table_db_manager = TableDBManager(
            self.db_temp_path_root,
            save_object['db_perm_path_root'],
            self.db_session_files
        )
        self.processed_table_manager = save_object['processed_table_manager']
        self.search_thread = SearchThread()
        self.search_preview_thread = FilePreviewThread()
        self.processing_thread = FileProcessingThread(self.table_db_manager)
        self.session_saved_flag = False
        self.n_parse_runs = save_object.get('n_parse_runs', 0)
        self.n_prunes = save_object.get('n_prunes', 0)
