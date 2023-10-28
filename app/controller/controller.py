from typing import TYPE_CHECKING, Callable
if TYPE_CHECKING:
    from uuid import UUID
    from model.model import Model
    from views.view import View
    from model.article_managers import (
        BaseData, Article, SuppFile, ProcessedTable
    )
    from views.page import (
        PageElements, SearchPageElements, ProcessedPageElements
    )
    from model.threading import SearchThread, FilePreviewThread, FileProcessingThread
    from views.list import DataListItem
    from PyQt5.QtWidgets import QListWidgetItem

from PyQt5.QtCore import Qt, QCoreApplication, QEventLoop
from PyQt5.QtWidgets import QFileDialog, QMessageBox

from utils.constants import PageIdentity, Mode


class Controller:
    def __init__(self, model: 'Model', view: 'View'):
        self.model = model
        self.view = view
        self._connect_sigs(
            view.search_elems,
            view.parsed_elems,
            view.pruned_elems,
            model.search_thread,
            model.search_preview_thread,
            model.processing_thread
        )

    def _set_state(self, state: 'Mode'):
        self.model.set_state(state)

    @property
    def curr_elems(self) -> 'PageElements':
        return self.view.active_elements

    @property
    def output_page(self) -> 'PageElements':
        mode: 'Mode' = self.model.state
        return {
            Mode.SEARCHING: self.view.search_elems,
            Mode.PROCESSING: self.view.parsed_elems,
            Mode.PRUNING: self.view.pruned_elems
        }.get(mode)

    def _connect_sigs(
            self,
            search: 'SearchPageElements',
            parsed: 'ProcessedPageElements',
            pruned: 'ProcessedPageElements',
            search_thread: 'SearchThread',
            search_preview_thread: 'FilePreviewThread',
            processing_thread: 'FileProcessingThread'
    ):
        signals_map = {
            self.view.save_action.triggered: self.save,
            self.view.load_action.triggered: self.load,

            # Search page
            search.article_ui_list.itemClicked: self.on_article_clicked,
            search.search_btn.clicked: self.search_articles,
            search.stop_search_btn.clicked: self.send_search_stop,
            search.proceed_btn.clicked: self.on_proceed,


            # Parsed page
            parsed.article_ui_list.itemClicked: self.on_article_clicked,
            parsed.filter_sig: self.filter_tables,
            parsed.prune_sig: self.prune_tables_and_columns,

            # Pruned page
            pruned.article_ui_list.itemClicked: self.on_article_clicked,
            pruned.filter_sig: self.filter_tables,
            pruned.prune_sig: self.prune_tables_and_columns,

            # Search threads
            search_thread.article_sig: self.display_article_in_list,
            search_thread.finished_sig: self.stop_search,
            search_preview_thread.prev_ready_sig: self.load_preview,

            # Processing thread
            processing_thread.article_sig: self.display_article_in_list,
            processing_thread.finished_sig: self.stop_processing
        }

        for signal, slot in signals_map.items():
            signal.connect(slot)

    def display_article_in_list(
        self,
        article: 'Article',
        progress: 'int',
        ids_list: 'list[tuple[str, UUID]]' = None
    ):
        if self.model.state == Mode.SEARCHING:
            article_data = self.model.create_article_data(article)
        elif self.model.state == Mode.PROCESSING:
            article_data = self.model.update_article(article, ids_list)

        self.view.display_article(self.output_page, article_data, progress)

    def search_articles(self):
        self._set_state(Mode.SEARCHING)
        self.model.reset_for_searching()
        self.view.tab_widget.setCurrentIndex(0)

        if self.model.search_thread.isRunning():
            QMessageBox.warning(
                self.view,
                "Search in Progress",
                "A search is already in progress. "
                "Please wait or stop the current search.")
            return

        query = self.output_page.query_field.text()
        if not query:
            return

        self.view.clear_list_and_observers(self.output_page.article_ui_list)
        self.view.show_searching_view()

        self.model.search_thread.prepare(query)
        self.model.search_thread.start()

    def stop_search(self, search_thread: 'SearchThread'):
        search_thread.quit()
        self.view.hide_searching_view()
        self._set_state(Mode.BROWSING)

    def send_search_stop(self):
        if self.model.search_thread.isRunning():
            self.model.search_thread.stop()
            self.model.search_thread.wait()

        while self.model.state == Mode.SEARCHING:
            QCoreApplication.processEvents(QEventLoop.AllEvents, 100)

    def stop_processing(self):
        if self.model.processing_thread.isRunning():
            self.model.processing_thread.stop()
            self.model.processing_thread.quit()
            self.model.processing_thread.wait()

        self.output_page.prog_bar.hide()
        self._set_state(Mode.BROWSING)

    def on_article_clicked(self, item: 'QListWidgetItem'):
        article_id = item.data(Qt.UserRole)
        article: 'Article' = self.model.bibliography.get_article(article_id)

        elements = self.curr_elems

        if elements.page_identity == PageIdentity.SEARCH:
            data_set = article.supp_files
        elif elements.page_identity == PageIdentity.PARSED:
            data_set = article.processed_tables
        elif elements.page_identity == PageIdentity.PRUNED:
            data_set = article.pruned_tables

        self.view.update_article_display(article, elements, data_set)

        for i in range(elements.data_ui_list.count()):
            list_item = elements.data_ui_list.item(i)
            widget: 'DataListItem' = elements.data_ui_list.itemWidget(
                list_item)
            if elements.page_identity == PageIdentity.SEARCH:
                widget.preview_requested.connect(self.request_suppfile_preview)
            else:
                widget.preview_requested.connect(self.preview_processed_table)

    def load_preview(
        self,
        data: 'BaseData',
        table_id=None,
        callback: 'Callable' = None
    ):

        # TODO idea: partially run the parser on the previewed file to get the
        # regions of the file that are tables, then highlight those regions in
        # the supp file previewer

        context = self.curr_elems.page_identity
        use_checkable_header = context != PageIdentity.SEARCH

        table: 'ProcessedTable' = self.model.processed_table_manager \
            .get_processed_table(table_id) if table_id else None

        cols = None
        if context == PageIdentity.PARSED:
            cols = table.checked_columns if table else None
        elif context == PageIdentity.PRUNED:
            cols = table.pruned_columns if table else None

        self.view.display_multisheet_table(
            data, use_checkable_header, table_id, callback, cols
        )

        self.view.stop_load_animation()

    # The original signal emits two arguments, but this slot only takes one
    # - why does this work? Granted, we don't need the context argument, but
    # still, it doesn't make sense why this doesn't throw an error.
    # TODO debug later
    def request_suppfile_preview(self, file_data: 'SuppFile'):
        self.view.start_load_animation()

        if self.model.search_preview_thread.isRunning():
            self.model.search_preview_thread.stop()
            self.model.search_preview_thread.quit()
            self.model.search_preview_thread.wait()

        # Need to manage the article:supp_file:table relationship such that
        # the metadata is available to the tables
        self.curr_elems.metadata_view.setHtml(file_data.metadata)

        self.model.search_preview_thread.prepare(file_data.url)
        self.model.search_preview_thread.start()

    def preview_processed_table(
        self, table: 'ProcessedTable', context: 'PageIdentity'
    ):
        table_data = {
            "sheet": self.model.table_db_manager.get_processed_table_data(
                table.id, context
            )
        }

        # Load the metadata; hack until I completely refactor the data model
        self.curr_elems.metadata_view.setHtml(table.supp_file.metadata)

        self.view.start_load_animation()
        self.load_preview(table_data, table.id, self._update_checked_columns)

    def _update_checked_columns(self, table_id, checked_columns):
        table: 'ProcessedTable' = self.model.processed_table_manager. \
            get_processed_table(table_id)

        if table:
            if self.curr_elems.page_identity == PageIdentity.PARSED:
                table.checked_columns = checked_columns
            elif self.curr_elems.page_identity == PageIdentity.PRUNED:
                table.pruned_columns = checked_columns

    def prune_tables_and_columns(self, context: 'PageIdentity'):
        self._set_state(Mode.PRUNING)

        for article in self.model.bibliography.articles.values():
            article.cascade_checked_state(context)

        self.view.set_active_tab(PageIdentity.PRUNED)
        self.model.prune_tables_and_columns(context)
        self.view.clear_page_lists(self.curr_elems)

        article: 'Article'
        for article in self.model.bibliography.get_selected_articles(context):
            self.view.display_article(self.output_page, article, 0)

        self._set_state(Mode.BROWSING)

    def filter_tables(self, context: 'PageIdentity'):
        query = self.curr_elems.query_filter_field.text()
        if not query:
            return
        self.model.filter_tables(query, context)

    def on_proceed(self):
        if self.model.search_thread.isRunning():
            reply = QMessageBox.question(
                self.view,
                "Search in Progress",
                "A search is still in progress. "
                "Do you want to stop the current search and proceed?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )

            if reply == QMessageBox.Yes:
                self.send_search_stop()
            else:
                return

        self.model.reset_for_processing()
        self.view.set_active_tab(PageIdentity.PARSED)

        if self.model.processing_thread.isRunning():
            QMessageBox.warning(
                self.view,
                "Processing in Progress",
                "A parsing run is already in progress. Please wait."
            )
            return

        self._set_state(Mode.PROCESSING)

        self.view.clear_list_and_observers(self.curr_elems.article_ui_list)
        self.curr_elems.prog_bar.setValue(0)
        self.curr_elems.prog_bar.show()

        for article in self.model.bibliography.articles.values():
            article.cascade_checked_state(PageIdentity.SEARCH)

        selected_articles = self.model.bibliography.get_selected_articles(
            PageIdentity.SEARCH
        )

        self.model.processing_thread.prepare(selected_articles)
        self.model.processing_thread.start()

    def save(self):
        # check if the app is in a state where it can be saved
        if self.model.state != Mode.BROWSING:
            QMessageBox.warning(
                self.view,
                "Cannot Save",
                "The application can only be saved in the browsing state."
            )
            return

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getSaveFileName(
            self.view,
            "Save As",
            f"{self.model.saves_path}/",  # Default location
            "All Files (*);;Text Files (*.txt)",
            options=options
        )

        if filename:
            self.model.save(filename)

    def load(self):
        if self.model.state != Mode.BROWSING:
            QMessageBox.warning(
                self.view,
                "Cannot Save",
                "The application can only be saved in the browsing state."
            )
            return

        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        filename, _ = QFileDialog.getOpenFileName(
            self.view,
            "Load File",
            "",
            "Pickle Files (*.pkl);;All Files (*)",
            options=options
        )

        if not filename:
            return

        self.model.load(filename, self.model.db_path, self.model.saves_path)

        # repopulate the GUI
        # clear all pages
        self.view.reset()

        # re-init
        # TODO this is a bad way to do this, should reset the model properly
        self.__init__(self.model, self.view)

        # Populate search page
        self._set_state(Mode.SEARCHING)
        selected_articles = self.model.bibliography.articles.values()
        for article in selected_articles:
            self.view.display_article(self.output_page, article, 0)

        self._set_state(Mode.BROWSING)

        # Populate parsed page
        if not self.model.n_parse_runs:
            return

        self._set_state(Mode.PROCESSING)
        selected_articles = self.model.bibliography.get_selected_articles(
            PageIdentity.SEARCH
        )

        for article in selected_articles:
            self.view.display_article(self.output_page, article, 0)

        self._set_state(Mode.BROWSING)

        # Populate pruned page
        if not self.model.n_prunes:
            return

        self._set_state(Mode.PRUNING)
        if self.model.n_prunes == 1:
            selected_articles = self.model.bibliography.get_selected_articles(
                PageIdentity.PARSED
            )

        else:
            selected_articles = self.model.bibliography.get_selected_articles(
                PageIdentity.PRUNED
            )

        for article in selected_articles:
            self.view.display_article(self.output_page, article, 0)

        self._set_state(Mode.BROWSING)
