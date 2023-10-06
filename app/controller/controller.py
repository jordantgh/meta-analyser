from PyQt5.QtCore import Qt, QCoreApplication, QEventLoop
from PyQt5.QtWidgets import QMessageBox
from enum import Enum, auto

class Mode(Enum):
    BROWSING = 0
    SEARCHING = auto()
    PROCESSING = auto()
    PRUNING = auto()

class Controller:
    def __init__(self, model, view):
        self.state = Mode.BROWSING
        self.model = model
        self.view = view
        self.search_page = self.view.search_components
        self.parsed_page = self.view.parsed_components
        self.pruned_page = self.view.pruned_components
        self.connect_sigs()

    def set_state(self, state):
        self.state = state
    
    @property
    def view_elem(self):
        return self.view.active_elements

    def connect_sigs(self):
        # Search page; article click
        self.search_page \
            .article_list.itemClicked \
            .connect(self.handle_article_click)

        # Search page; buttons
        self.search_page \
            .search_btn.clicked \
            .connect(self.search_articles)
        self.search_page \
            .stop_search_btn.clicked \
            .connect(self.send_stop)
        self.search_page \
            .proceed_btn.clicked \
            .connect(self.on_proceed)

        # Search thread
        self.model \
            .search_thread.article_sig \
            .connect(self.add_article)
        self.model \
            .search_thread.finished_sig \
            .connect(self.stop_search)

        # Processing thread
        self.model \
            .processing_thread.article_sig \
            .connect(self.on_article_processed)
        self.model \
            .processing_thread.finished_sig \
            .connect(self.on_processing_finished)
        self.model \
            .preview_thread.prev_ready_sig \
            .connect(self.load_preview)

        # Parsed tables page; article click
        self.parsed_page \
            .article_list.itemClicked \
            .connect(self.handle_processed_article_click)

        # Parsed tables page; buttons
        self.parsed_page \
            .filter_btn.clicked \
            .connect(self.filter_tables)
        self.parsed_page \
            .prune_btn.clicked \
            .connect(self.prune_tables_and_columns)

        # Pruned tables page; article click
        self.pruned_page \
            .article_list.itemClicked \
            .connect(self.handle_pruned_article_click)

        # Pruned tables page; buttons
        self.pruned_page \
            .filter_btn.clicked \
            .connect(self.filter_tables)
        self.pruned_page \
            .prune_btn.clicked \
            .connect(self.prune_tables_and_columns)

    def on_article_discovered(self, article_json, progress):
        article_data = self.model.create_article_data(article_json)
        self.view.display_article(self.search_page, article_data, progress)

    def on_article_processed(self, article, ids_list, progress):
        article_data = self.model.update_article(article, ids_list)
        self.view.display_article(self.parsed_page, article_data, progress)

    def search_for_articles(self):
        self.set_state(Mode.SEARCHING)
        self.model.reset_for_searching()
        self.view.tab_widget.setCurrentIndex(0)
        if self.model.search_thread.isRunning():
            QMessageBox.warning(
                self.view,
                "Search in Progress",
                "A search is already in progress. " 
                "Please wait or stop the current search.")
            return

        self.model.search_thread.should_stop = False
        query = self.search_page.query_field.text()
        if not query:
            return
        
        # TODO these are low level concerns that should be handled by the view
        self.search_page.article_list.clear()
        self.search_page.previews.hide()
        self.search_page.prog_bar.setValue(0)
        self.search_page.prog_bar.show()
        self.search_page.search_status.setText("Searching...")
        self.search_page.stop_search_btn.show()
        self.search_page.stop_search_btn.setEnabled(True)
        
        self.model.search_thread.query = query
        self.model.search_thread.start()

    def stop_search(self, search_thread):
        search_thread.quit()
        
        # TODO these are low level concerns that should be handled by the view
        self.search_page.search_status.setText("Stopping search...")
        self.search_page.prog_bar.hide()        
        self.search_page.search_status.setText("Search stopped.")
        self.search_page.stop_search_btn.hide()
        self.search_page.stop_search_btn.setEnabled(False)
        
        self.set_state(Mode.BROWSING)
        
    def send_search_stop(self):
        if self.model.search_thread.isRunning():            
            self.model.search_thread.stop()
            self.model.search_thread.wait()

        while self.state == Mode.SEARCHING:
            QCoreApplication.processEvents(QEventLoop.AllEvents, 100)

    def on_processing_finished(self):
        if self.model.processing_thread.isRunning():
            self.model.processing_thread.stop()
            self.model.processing_thread.quit()
            self.model.processing_thread.wait()

        self.parsed_page.prog_bar.hide()
        self.set_state(Mode.BROWSING)

    # TODO these three can be combined and list item type passed in or derived
    # from view_elem.__class__.__name__
    def handle_article_click(self, item):
        self.view_elem.previews.hide()
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'supp_files',
            self.view.suppfilelistitem_factory)
        
        # TODO do we need to pass in the factory if the relevant class can be
        # derived from list item type within view.update_article_display?
        # ... could just do away with the factories altogether

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_supp_file)

    def handle_processed_article_click(self, item):
        self.view_elem.previews.hide()
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'processed_tables',
            self.view.processedtablelistitem_factory)

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_processed_table)

    def handle_pruned_article_click(self, item):
        self.view_elem.previews.hide()
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)
        
        self.view.update_article_display(
            article,
            'pruned_article_tables',
            self.view.processedtablelistitem_factory)

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_pruned_table)

    def load_preview(self, data, table_id=None, callback=None):
        use_checkable_header = self.view_elem \
            .__class__.__name__ == 'ProcessedPageElements'

        for i in reversed(range(self.view_elem.previews_layout.count())):
            widget = self.view_elem.previews_layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        processed_table = self.model.processed_table_manager \
            .get_processed_table(table_id) if table_id else None

        checked_columns = processed_table \
            .checked_columns if processed_table else None

        self.view.display_multisheet_table(
            data, use_checkable_header, table_id, callback, checked_columns)
        self.view.stop_load_animation()
        self.view_elem.previews.show()

    def preview_supp_file(self, file_id):
        file_data = self.model.file_manager.get_file(file_id)
        self.view.start_load_animation()
        
        if self.model.preview_thread.isRunning():
            self.model.preview_thread.quit()
            self.model.preview_thread.wait()
            
        self.model.preview_thread.file_url = file_data.url
        self.model.preview_thread.start()

    def preview_processed_table(self, table_id):
        table_data = {
            "sheet": self.model.table_db_manager.get_processed_table_data(
                table_id)}
        
        self.view.start_load_animation()
        self.load_preview(table_data, table_id, self.update_checked_columns)

    def preview_pruned_table(self, table_id):
        table_data = {
            "sheet": self.model.table_db_manager.get_post_pruning_table_data(
                table_id)}
        
        self.view.start_load_animation()
        self.load_preview(table_data, table_id, self.update_checked_columns)

    def update_checked_columns(self, table_id, checked_columns):
        processed_table = self.model.processed_table_manager \
            .get_processed_table(table_id)

        if processed_table:
            processed_table.checked_columns = checked_columns

    def prune_tables_and_columns(self):
        self.view.tab_widget.setCurrentIndex(2)
        self.model.prune_tables_and_columns()
        self.view.clear_article_list_and_files_view()
        
        # display all selected articles in the pruned page
        # TODO need to make it so articles w/ no tables dont appear here
        for article in self.model.bibliography.get_selected_articles():
            self.view.display_article(self.pruned_page, article, 0)
    def filter_tables(self):
        query = self.view_elem.query_filter_field.text()
        if not query:
            return
        self.model.filter_tables(query)

    def on_proceed(self):
        if self.model.search_thread.isRunning():
            reply = QMessageBox.question(
                self.view,
                "Search in Progress",
                "A search is still in progress. "
                "Do you want to stop the current search and proceed?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No)
            
            if reply == QMessageBox.Yes:
                self.send_search_stop()
            else:
              return

        self.model.reset_for_processing()
        self.view.tab_widget.setCurrentIndex(1)
        
        if self.model.processing_thread.isRunning():
            QMessageBox.warning(
                self.view,
                "Processing in Progress",
                "A parsing run is already in progress. Please wait.")
            return

        self.set_state(Mode.PROCESSING)
        self.model.processing_thread.should_stop = False
        
        # TODO these are low level concerns that should be handled by the view
        self.view_elem.article_list.clear()
        self.view_elem.previews.hide()
        self.view_elem.prog_bar.setValue(0)
        self.view_elem.prog_bar.show()
        
        selected_articles = self.model.bibliography.get_selected_articles()
        self.model.processing_thread.selected_articles = selected_articles
        self.model.processing_thread.start()
