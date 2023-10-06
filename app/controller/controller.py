from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox


class Controller:
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.search_page = self.view.search_components
        self.parsed_page = self.view.parsed_components
        self.pruned_page = self.view.pruned_components
        self.connect_sigs()

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


    def add_article(self, article_json, progress):
        article_data = self.model.create_article_data(article_json)
        self.view.display_article(article_data, progress)

    def on_article_processed(self, article, ids_list, progress):
        article_data = self.model.update_article(article, ids_list)
        self.view.display_article(article_data, progress)

    def search_articles(self):
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
        query = self.view_elem.query_field.text()
        if not query:
            return
        
        # these are low level view details that should be handled by the view
        self.view_elem.article_list.clear()
        self.view_elem.previews.hide()
        self.view_elem.prog_bar.setValue(0)
        self.view_elem.prog_bar.show()
        self.view_elem.search_status.setText("Searching...")
        
        self.model.search_thread.query = query
        self.model.search_thread.start()

        self.view_elem.stop_search_btn.show()
        self.view_elem.stop_search_btn.setEnabled(True)

    def stop_search(self):
        if self.model.search_thread.isRunning():
            self.model.search_thread.stop()
            self.model.search_thread.quit()
            self.view_elem.search_status.setText("Stopping search...")
            self.model.search_thread.wait()
        self.view_elem.prog_bar.hide()
        self.view_elem.search_status.setText("Search stopped.")
        
        self.view_elem.stop_search_btn.hide()
        self.view_elem.stop_search_btn.setEnabled(False)

    def on_search_finished(self):
        self.view_elem.prog_bar.hide()        
        self.view_elem.search_status.clear()
        self.view_elem.stop_search_btn.hide()
        self.view_elem.stop_search_btn.setEnabled(False)

    def on_processing_finished(self):
        self.view_elem.prog_bar.hide()  

    def handle_article_click(self, item):
        self.view_elem.previews.hide()
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'supp_files',
            self.view.suppfilelistitem_factory)

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

    def refresh_view(self):
        self.view.clear_article_list_and_files_view()
        
        for article in self.model.bibliography.get_selected_articles():
            self.view.display_article(article, 0)

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
        self.refresh_view()

    def filter_tables(self):
        query = self.view_elem.query_filter_field.text()
        if not query:
            return
        self.model.filter_tables(query)

    def on_proceed(self):
        self.model.reset_for_processing()
        self.view.tab_widget.setCurrentIndex(1)
        
        if self.model.processing_thread.isRunning():
            QMessageBox.warning(
                self.view,
                "Processing in Progress",
                "A parsing run is already in progress. Please wait.")
            return
        self.model.processing_thread.should_stop = False
        self.view_elem.article_list.clear()
        self.view_elem.previews.hide()
        self.view_elem.prog_bar.setValue(0)
        self.view_elem.prog_bar.show()
        
        selected_articles = self.model.bibliography.get_selected_articles()
        self.model.processing_thread.selected_articles = selected_articles
        self.model.processing_thread.start()
