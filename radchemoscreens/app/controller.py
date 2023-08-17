import json
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QMessageBox
from model import SearchThread, FilePreviewThread, Bibliography

class CRISPRController:
    def __init__(self, model, view):
        self.model = model
        self.view = view
        self.view.init_load_animation()
        self.bibliography = Bibliography()
        self.connect_sigs()

    def connect_sigs(self):
        self.view.search_btn.clicked.connect(self.search_papers)
        self.view.paper_list.itemClicked.connect(self.handle_paper_click)

        self.search_thread = SearchThread()
        self.search_thread.article_sig.connect(self.add_paper)
        self.search_thread.finished_sig.connect(self.on_search_finished)
        self.preview_thread = FilePreviewThread("")
        self.preview_thread.prev_ready_sig.connect(self.load_preview)
        
        self.view.stop_search_btn.clicked.connect(self.stop_search)
        self.view.proceed_btn.clicked.connect(self.on_proceed)

    def add_paper(self, article, progress):
        paper_data = self.model.create_paper_data(article)
        self.bibliography.add_paper(paper_data)
        self.view.display_paper_in_ui(paper_data, progress)

    def search_papers(self):
        if self.search_thread.isRunning():
            QMessageBox.warning(self.view, "Search in Progress", "A search is already in progress. Please wait or stop the current search.")
            return
        self.search_thread.should_stop = False
        query = self.view.query_field.text()
        if not query: return
        self.view.paper_list.clear()
        self.view.prog_bar.setValue(0)
        self.view.prog_bar.show()
        self.view.search_status.setText("Searching...")
        self.search_thread.query = query
        self.search_thread.start()
        
        self.view.stop_search_btn.setEnabled(True)

    def stop_search(self):
        if self.search_thread.isRunning():
            self.search_thread.stop()
            self.search_thread.quit()
            self.search_thread.wait()
            self.view.search_status.setText("Stopping search...")
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
        
        for i in range(self.view.supp_files.count()):
            list_item = self.view.supp_files.item(i)
            widget = self.view.supp_files.itemWidget(list_item)
            widget.preview_requested.connect(self.preview_supp_file)

    def load_preview(self, fname):
        data = self.model.extract_data_from_file(fname)
        self.view.stop_load_animation()

        for i in reversed(range(self.view.previews_layout.count())):
            self.view.previews_layout.itemAt(i).widget().setParent(None)

        self.determine_data_type_and_display(data)
        self.view.previews.show()

    def preview_supp_file(self, file_url):
        self.view.start_load_animation()
        if self.preview_thread.isRunning():
            self.preview_thread.quit()
            self.preview_thread.wait()
        self.preview_thread.file_url = file_url
        self.preview_thread.start()

    def determine_data_type_and_display(self, data):
        if isinstance(data, dict):
            self.view.display_data_in_tabs(data)
        else:
            self.view.display_data_in_table(data)

    def on_proceed(self):
        selected_data = self.extract_selected_papers_and_files()
        with open('selected_papers_and_files.json', 'w') as f:
            json.dump(selected_data, f, indent=4)
        QMessageBox.information(self.view, "Saved", "Selected papers and files have been saved to selected_papers_and_files.json")

    def extract_selected_papers_and_files(self):
        selected_papers = {}
        for paper in self.bibliography.get_selected_papers():
            selected_files = [f.url for f in paper.files if f.checked]
            if selected_files:
                selected_papers[paper.pmc_id] = {
                    "title": paper.title,
                    "files": selected_files
                }
        return selected_papers