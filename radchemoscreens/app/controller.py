import json
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox


class Controller:
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.connect_sigs()

    def connect_sigs(self):
        self.view.search_btn.clicked.connect(self.search_articles)
        self.view.article_list.itemClicked.connect(self.handle_article_click)
        self.view.stop_search_btn.clicked.connect(self.stop_search)
        self.view.proceed_btn.clicked.connect(self.on_proceed)

        self.model.search_thread.article_sig.connect(self.add_article)
        self.model.search_thread.finished_sig.connect(self.on_search_finished)

        self.model.processing_thread.article_sig.connect(self.on_article_processed)
        self.model.processing_thread.finished_sig.connect(self.on_search_finished)

        self.model.preview_thread.prev_ready_sig.connect(self.load_preview)

    def add_article(self, article, progress):
        article_data = self.model.create_article_data(article)
        self.model.bibliography.add_article(article_data)
        self.view.display_article(article_data, progress)
        
    def on_article_processed(self, article, ids_list, progress):
        article_data = self.model.update_article(article, ids_list)
        self.view.display_article(article_data, progress)

    def search_articles(self):
        self.processing_mode = False
        if self.model.search_thread.isRunning():
            QMessageBox.warning(self.view, "Search in Progress", "A search is already in progress. Please wait or stop the current search.")
            return
        self.model.search_thread.should_stop = False
        query = self.view.query_field.text()
        if not query: return
        self.view.article_list.clear()
        self.view.previews.hide()
        self.view.prog_bar.setValue(0)
        self.view.prog_bar.show()
        self.view.search_status.setText("Searching...")
        self.model.search_thread.query = query
        self.model.search_thread.start()
        
        self.view.stop_search_btn.setEnabled(True)

    def stop_search(self):
        if self.model.search_thread.isRunning():
            self.model.search_thread.stop()
            self.model.search_thread.quit()
            self.view.search_status.setText("Stopping search...")
            self.model.search_thread.wait()
        self.view.prog_bar.hide()
        self.view.search_status.setText("Search stopped.")
        self.view.stop_search_btn.setEnabled(False)

    def on_search_finished(self):
        print("Finished.")
        self.view.prog_bar.hide()
        self.view.search_status.clear()
        self.view.stop_search_btn.setEnabled(False)

    def handle_article_click(self, item):
        self.view.previews.hide()
        article_id = item.data(Qt.UserRole)
        article = self.model.bibliography.get_article(article_id)
        
        if self.processing_mode:
            self.view.update_article_display_processed(article)
        else:
            self.view.update_article_display(article)

        for i in range(self.view.supp_files_view.count()):
            list_item = self.view.supp_files_view.item(i)
            widget = self.view.supp_files_view.itemWidget(list_item)
            if self.processing_mode:
                widget.preview_requested.connect(self.preview_processed_table)
            else:
                widget.preview_requested.connect(self.preview_supp_file)

    def load_preview(self, data):

        for i in reversed(range(self.view.previews_layout.count())):
            widget = self.view.previews_layout.itemAt(i).widget()
            if widget: widget.deleteLater()

        self.view.display_multisheet_table(data)
        self.view.stop_load_animation()
        self.view.previews.show()

    def preview_supp_file(self, file_id):
        file_data = self.model.file_manager.get_file(file_id)
        self.view.start_load_animation()
        if self.model.preview_thread.isRunning():
            self.model.preview_thread.quit()
            self.model.preview_thread.wait()
        self.model.preview_thread.file_url = file_data.url
        self.model.preview_thread.start()
        
    def preview_processed_table(self, table_id):
        print(f"previewing processed table {table_id}")
        table_data = {"sheet": self.model.table_db_manager.get_table(table_id)}
        print(f"preview of table (of type {type(table_data)}): {table_data}")
        self.view.start_load_animation()
        self.load_preview(table_data)
        
    def on_proceed(self):
        self.processing_mode = True
        if self.model.processing_thread.isRunning():
            QMessageBox.warning(self.view, "Processing in Progress", "A parsing run is already in progress. Please wait.")
            return
        self.model.processing_thread.should_stop = False
        query = self.view.query_field.text()
        if not query: return
        self.view.article_list.clear()
        self.view.previews.hide()
        self.view.prog_bar.setValue(0)
        self.view.prog_bar.show()
        selected_articles = self.model.bibliography.get_selected_articles()
        self.model.processing_thread.selected_articles = selected_articles
        self.model.processing_thread.start()

        # minimal debug code:
        all_selected_data = self.model.extract_selected_articles_and_files()
        with open('selected_articles_and_files.json', 'w') as f:
            json.dump(all_selected_data, f, indent=4)
        QMessageBox.information(self.view, "Saved", "Selected articles and files have been saved to selected_articles_and_files.json")