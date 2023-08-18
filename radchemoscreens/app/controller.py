import json
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox

class Controller:
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.connect_sigs()

    def connect_sigs(self):
        self.view.search_btn.clicked.connect(self.search_papers)
        self.view.paper_list.itemClicked.connect(self.handle_paper_click)
        self.view.stop_search_btn.clicked.connect(self.stop_search)
        self.view.proceed_btn.clicked.connect(self.on_proceed)

        self.model.search_thread.article_sig.connect(self.add_paper)
        self.model.search_thread.finished_sig.connect(self.on_search_finished)
        self.model.preview_thread.prev_ready_sig.connect(self.load_preview)
        
    def add_paper(self, article, progress):
        paper_data = self.model.create_paper_data(article)
        self.model.bibliography.add_paper(paper_data)
        self.view.display_paper(paper_data, progress)

    def search_papers(self):
        if self.model.search_thread.isRunning():
            QMessageBox.warning(self.view, "Search in Progress", "A search is already in progress. Please wait or stop the current search.")
            return
        self.model.search_thread.should_stop = False
        query = self.view.query_field.text()
        if not query: return
        self.view.paper_list.clear()
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
        self.view.prog_bar.hide()
        self.view.search_status.clear()
        self.view.stop_search_btn.setEnabled(False)

    def handle_paper_click(self, item):
        paper = item.data(Qt.UserRole)
        self.view.update_paper_display(paper)
        
        for i in range(self.view.supp_files_view.count()):
            list_item = self.view.supp_files_view.item(i)
            widget = self.view.supp_files_view.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_supp_file)

    def load_preview(self, data):
        self.view.stop_load_animation()

        for i in reversed(range(self.view.previews_layout.count())):
            widget = self.view.previews_layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        self.determine_data_type_and_display(data)
        self.view.previews.show()

    def preview_supp_file(self, file_url):
        self.view.start_load_animation()
        if self.model.preview_thread.isRunning():
            self.model.preview_thread.quit()
            self.model.preview_thread.wait()
        self.model.preview_thread.file_url = file_url
        self.model.preview_thread.start()

    def determine_data_type_and_display(self, data):
        if isinstance(data, dict):
            self.view.display_multisheet_table(data)
        else:
            self.view.display_table(data)

    def on_proceed(self):
        selected_data = self.model.extract_selected_papers_and_files()
        with open('selected_papers_and_files.json', 'w') as f:
            json.dump(selected_data, f, indent=4)
        QMessageBox.information(self.view, "Saved", "Selected papers and files have been saved to selected_papers_and_files.json")