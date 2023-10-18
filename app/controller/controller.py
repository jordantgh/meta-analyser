from PyQt5.QtCore import Qt, QCoreApplication, QEventLoop
from PyQt5.QtWidgets import QFileDialog, QMessageBox


class Controller:
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.search_page = self.view.search_components
        self.parsed_page = self.view.parsed_components
        self.pruned_page = self.view.pruned_components
        self.connect_sigs()

    def set_state(self, state):
        self.model.set_state(state)

    @property
    def view_elem(self):
        return self.view.active_elements

    def connect_sigs(self):
        # Menu bar
        self.view.save_action.triggered.connect(self.save)
        self.view.load_action.triggered.connect(self.load)

        # Search page; article click
        self.search_page \
            .article_list.itemClicked \
            .connect(self.handle_article_click)
        # Search page; buttons
        self.search_page \
            .search_btn.clicked \
            .connect(self.search_for_articles)
        self.search_page \
            .stop_search_btn.clicked \
            .connect(self.send_search_stop)
        self.search_page \
            .proceed_btn.clicked \
            .connect(self.on_proceed)
        # Search thread
        self.model \
            .search_thread.article_sig \
            .connect(self.on_article_discovered)
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
        self.view.display_article(
            self.search_page, 'search', article_data, progress)

    def on_article_processed(self, article, ids_list, progress):
        article_data = self.model.update_article(article, ids_list)
        self.view.display_article(
            self.parsed_page, 'parsed', article_data, progress)

    def search_for_articles(self):
        self.set_state(self.model.Mode.SEARCHING)
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
        self.view.clear_list_and_observers(self.search_page.article_list)
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

        self.set_state(self.model.Mode.BROWSING)

    def send_search_stop(self):
        if self.model.search_thread.isRunning():
            self.model.search_thread.stop()
            self.model.search_thread.wait()

        while self.model.state == self.model.Mode.SEARCHING:
            QCoreApplication.processEvents(QEventLoop.AllEvents, 100)

    def on_processing_finished(self):
        if self.model.processing_thread.isRunning():
            self.model.processing_thread.stop()
            self.model.processing_thread.quit()
            self.model.processing_thread.wait()

        self.parsed_page.prog_bar.hide()
        self.set_state(self.model.Mode.BROWSING)

    # TODO these three can be combined and list item type passed in or derived
    # from view_elem.__class__.__name__
    def handle_article_click(self, item):
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'supp_files',
            self.view.suppfilelistitem_factory,
            'search')

        # TODO do we need to pass in the factory if the relevant class can be
        # derived from list item type within view.update_article_display?
        # ... could just do away with the factories altogether

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_supp_file)

    def handle_processed_article_click(self, item):
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'processed_tables',
            self.view.processedtablelistitem_factory,
            'parsed')

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_processed_table)

    def handle_pruned_article_click(self, item):
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)

        self.view.update_article_display(
            article,
            'pruned_article_tables',
            self.view.processedtablelistitem_factory,
            'pruned')

        for i in range(self.view_elem.supp_files_view.count()):
            list_item = self.view_elem.supp_files_view.item(i)
            widget = self.view_elem.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_pruned_table)

    def load_preview(self, data, table_id=None, callback=None, context=None):
        use_checkable_header = self.view_elem \
            .__class__.__name__ == 'ProcessedPageElements'

        processed_table = self.model.processed_table_manager \
            .get_processed_table(table_id) if table_id else None

        checked_columns = None
        if context == 'parsed':
            checked_columns = processed_table \
                .checked_columns if processed_table else None
        elif context == 'pruned':
            checked_columns = processed_table \
                .pruned_columns if processed_table else None

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
                table_id).head(100)}

        self.view.start_load_animation()
        self.load_preview(table_data, table_id,
                          self.update_checked_columns, 'parsed')

    def preview_pruned_table(self, table_id):
        table_data = {
            "sheet": self.model.table_db_manager.get_post_pruning_table_data(
                table_id).head(100)}

        self.view.start_load_animation()
        self.load_preview(table_data, table_id,
                          self.update_pruned_columns, 'pruned')

    def update_checked_columns(self, table_id, checked_columns):
        processed_table = self.model.processed_table_manager \
            .get_processed_table(table_id)

        if processed_table:
            processed_table.checked_columns = checked_columns

    def update_pruned_columns(self, table_id, checked_columns):
        processed_table = self.model.processed_table_manager \
            .get_processed_table(table_id)

        if processed_table:
            processed_table.pruned_columns = checked_columns

    def prune_tables_and_columns(self):
        self.set_state(self.model.Mode.PRUNING)
        # get current page from view_elem (class name doesnt work!)
        current_page = self.view.tab_widget.currentIndex()
        if current_page == 1:
            context = 'parsed'
        elif current_page == 2:
            context = 'pruned'

        for article in self.model.bibliography.articles.values():
            article.cascade_checked_state(context)

        self.view.tab_widget.setCurrentIndex(2)
        self.model.prune_tables_and_columns(context)

        self.view.clear_article_list_and_files_view()

        # display all selected articles in the pruned page
        for article in self.model.bibliography.get_selected_articles(context):
            self.view.display_article(self.pruned_page, 'pruned', article, 0)

        self.set_state(self.model.Mode.BROWSING)
        
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

        self.set_state(self.model.Mode.PROCESSING)
        self.model.processing_thread.should_stop = False

        # TODO these are low level concerns that should be handled by the view
        self.view.clear_list_and_observers(self.view_elem.article_list)
        self.view_elem.prog_bar.setValue(0)
        self.view_elem.prog_bar.show()

        for article in self.model.bibliography.articles.values():
            article.cascade_checked_state('search')

        selected_articles = self.model.bibliography.get_selected_articles(
            'search')
        self.model.processing_thread.selected_articles = selected_articles
        self.model.processing_thread.start()

    def save(self):
        # check if the app is in a state where it can be saved
        if self.model.state != self.model.Mode.BROWSING:
            QMessageBox.warning(
                self.view,
                "Cannot Save",
                "The application can only be saved in the browsing state.")
            return

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getSaveFileName(
            self.view, "Save As", "", "All Files (*);;Text Files (*.txt)", options=options)

        if filename:
            self.model.save(filename)
            
    def load(self):        
        if self.model.state != self.model.Mode.BROWSING:
            QMessageBox.warning(
                self.view,
                "Cannot Save",
                "The application can only be saved in the browsing state.")
            return

        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        filename, _ = QFileDialog.getOpenFileName(
            self.view, "Load File", "", "Pickle Files (*.pkl);;All Files (*)", options=options
        )

        if not filename:
            return

        self.model.load(filename)
        
        # repopulate the GUI
        ## clear all pages
        for page in [self.search_page, self.parsed_page, self.pruned_page]:
            self.view.clear_list_and_observers(page.article_list)
            page.supp_files_view.clear()
            page.title_abstract_disp.clear()
            page.previews.clear()
            
        # re-init 
        self.__init__(self.model, self.view)
        
        ## populate search page
        selected_articles = self.model.bibliography.articles.values()
        
        for article in selected_articles:
            self.view.display_article(self.search_page, 'search', article, 0)
        
        ## populate parsed page
        if not self.model.ever_parsed:
            return
        
        selected_articles = self.model.bibliography \
            .get_selected_articles('search')
            
        for article in selected_articles:
            self.view.display_article(self.parsed_page, 'parsed', article, 0)
            
        ## populate pruned page
        
        # BUG: the current implementation wont handle it if there was >1
        # prune run in the loaded model. If pruned >1 times, we should get
        # selected articles from the pruned page and not the parsed page
        print(self.model.ever_pruned)
        if not self.model.ever_pruned:
            return
        
        if self.model.ever_pruned == 1:
            selected_articles = self.model.bibliography \
                .get_selected_articles('parsed')
        else:
            selected_articles = self.model.bibliography \
                .get_selected_articles('pruned')

        for article in selected_articles:
            self.view.display_article(self.pruned_page, 'pruned', article, 0)