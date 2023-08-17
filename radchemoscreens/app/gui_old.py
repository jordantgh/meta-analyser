import os
import sys
import logging
import requests
import json
import pandas as pd

from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QPushButton, QListWidget, QListWidgetItem, QProgressBar, QTextEdit, QCheckBox, QPushButton, QTableWidget, QTableWidgetItem, QSizePolicy, QTabWidget, QMessageBox

from search_for_papers import query_pmc

class SearchThread(QThread):
    article_sig = pyqtSignal(dict, int)
    finished_sig = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.query = ""
        self.should_stop = False

    def stop(self):
        self.should_stop = True

    def run(self):
        try:
            query_pmc(self.query, callback=self.article_sig.emit, thread = self)
            self.finished_sig.emit()
        except Exception as e:
            logging.error(f"Unhandled exception in SearchThread: {e}")

class FilePreviewThread(QThread):
    prev_ready_sig = pyqtSignal(str)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url

    def run(self):
        fname = self.download_file(self.file_url)
        self.prev_ready_sig.emit(fname)
        
    def download_file(self, url):
        local_fname = url.split('/')[-1]
        if not local_fname.endswith(('.txt', '.csv', '.xlsx')):
            return
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_fname, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return local_fname

class SupplementaryFile:
    def __init__(self, url):
        self.url = url
        self.checked = True

class Paper:
    def __init__(self, title, abstract, pmc_id, files=[]):
        self.title = title
        self.abstract = abstract
        self.pmc_id = pmc_id
        self.files = [SupplementaryFile(url) for url in files]
        self.checked = True

class Bibliography:
    def __init__(self):
        self.papers = []

    def add_paper(self, paper):
        self.papers.append(paper)
        
    def get_selected_papers(self):
        return [p for p in self.papers if p.checked]

class UIListItem(QWidget):
    def __init__(self, title, data=None):
        super().__init__()
        self.data = data
        self.checkbox = QCheckBox()
        self.checkbox.setChecked(self.data.checked)
        self.checkbox.toggled.connect(self.checkbox_toggled)
        
        label = QLabel(title)
        label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        layout = QHBoxLayout(self)
        layout.addWidget(self.checkbox)
        layout.addWidget(label)
        layout.addStretch(1)
        self.setLayout(layout)

    def checkbox_toggled(self):
        self.data.checked = self.checkbox.isChecked()

    # override default mousePressEvent behaviour
    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)

class SuppFileItem(UIListItem):
    def __init__(self, file_url, main_window, file_instance):
        self.file_url = file_url
        self.main_window = main_window
        disp_name = self.get_disp_name(file_url.split('/')[-1])

        super().__init__(disp_name, file_instance)

        self.preview_btn = QPushButton("Preview", self)
        self.preview_btn.clicked.connect(self.preview_file)
        self.layout().addWidget(self.preview_btn)

    def get_disp_name(self, text):
        font_metrics = QFontMetrics(self.main_window.font())
        available_width = self.main_window.supp_files.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)

    def preview_file(self):
        self.main_window.preview_supp_file(self.file_url)

class CRISPRApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        self.init_uis()
        self.setup_layouts()
        self.connect_sigs()
        self.bibliography = Bibliography()
        self.init_load_animation()

    def init_uis(self):
        # Central widget
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        self.stop_search_btn = QPushButton("Stop Search", self)
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", self)

        # Paper list elements
        self.query_field = QLineEdit(self)
        self.search_status = QLabel("")
        self.search_btn = QPushButton("Search", self)
        self.paper_list = QListWidget(self)
        self.prog_bar = QProgressBar(self)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()

        # Paper details elements
        self.title_disp = QTextEdit(self)
        self.title_disp.setPlaceholderText("Title will be shown here")
        self.abstract_disp = QTextEdit(self)
        self.abstract_disp.setPlaceholderText("Abstract will be shown here")
        self.supp_files = QListWidget(self)
        self.load_label = QLabel(self)
        self.load_label.setAlignment(Qt.AlignCenter)
        self.previews = QWidget(self)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)

    def setup_layouts(self):
        # Paper list layout
        pane_0 = QVBoxLayout()
        pane_0.addWidget(QLabel("Enter Query:"))
        pane_0.addWidget(self.query_field)
        pane_0.addWidget(self.prog_bar)
        pane_0.addWidget(self.search_btn)
        pane_0.addWidget(self.stop_search_btn)
        pane_0.addWidget(self.paper_list)

        # Paper details layout
        pane_1 = QVBoxLayout()
        pane_1.addWidget(QLabel("Title:"))
        pane_1.addWidget(self.title_disp)
        pane_1.addWidget(QLabel("Abstract:"))
        pane_1.addWidget(self.abstract_disp)
        pane_1.addWidget(QLabel("Supplementary Files:"))
        pane_1.addWidget(self.supp_files)
        pane_1.addWidget(self.load_label)

        # Main layout
        main_pane = QHBoxLayout(self.centralWidget())
        main_pane.addLayout(pane_0)
        main_pane.addLayout(pane_1)
        main_pane.addWidget(self.previews)

    def connect_sigs(self):
        self.search_btn.clicked.connect(self.search_papers)
        self.paper_list.itemClicked.connect(self.handle_paper_click)

        # Thread signals
        self.search_thread = SearchThread()
        self.search_thread.article_sig.connect(self.add_paper)
        self.search_thread.finished_sig.connect(self.on_search_finished)

        self.preview_thread = FilePreviewThread("")
        self.preview_thread.prev_ready_sig.connect(self.load_preview)
        
        self.stop_search_btn.clicked.connect(self.stop_search)
        self.proceed_btn.clicked.connect(self.on_proceed)

    def init_load_animation(self):
        self.load_timer = QTimer(self)
        self.load_dots = 0
        self.load_timer.timeout.connect(self.update_load_text)

    def start_load_animation(self):
        self.load_dots = 0
        self.load_label.setText("Loading.")
        self.load_timer.start(500)  # ms

    def stop_load_animation(self):
        self.load_timer.stop()
        self.load_label.clear()

    def update_load_text(self):
        self.load_dots = (self.load_dots + 1) % 4
        self.load_label.setText("Loading" + "." * self.load_dots)

    def add_paper(self, article, progress):
        paper_data = self.create_paper_data(article)
        self.bibliography.add_paper(paper_data)
        self.display_paper_in_ui(paper_data, progress)

    def create_paper_data(self, article):
        return Paper(
            title=article["Title"],
            abstract=article["Abstract"],
            pmc_id=article["PMCID"],
            files=article["SupplementaryFiles"]
        )

    def display_paper_in_ui(self, paper_data, progress):
        item = QListWidgetItem()
        paper_name = UIListItem(paper_data.title, paper_data)
        item.setSizeHint(paper_name.sizeHint())
        item.setData(Qt.UserRole, paper_data)
        self.paper_list.addItem(item)
        self.paper_list.setItemWidget(item, paper_name)
        self.prog_bar.setValue(progress + 1)

    def search_papers(self):
        if self.search_thread.isRunning():
            QMessageBox.warning(self, "Search in Progress", "A search is already in progress. Please wait or stop the current search.")
            return
        self.search_thread.should_stop = False
        query = self.query_field.text()
        if not query: return
        self.paper_list.clear()
        self.prog_bar.setValue(0)
        self.prog_bar.show()
        self.search_status.setText("Searching...")
        self.search_thread.query = query
        self.search_thread.start()
        
        self.stop_search_btn.setEnabled(True)

    def stop_search(self):
        if self.search_thread.isRunning():
            self.search_thread.stop()
            self.search_thread.quit()
            self.search_thread.wait()
            self.search_status.setText("Stopping search...")
        self.prog_bar.hide()
        self.search_status.setText("Search stopped.")
        self.stop_search_btn.setEnabled(False)

    def on_search_finished(self):
        self.prog_bar.hide()
        self.search_status.clear()
        self.stop_search_btn.setEnabled(False)

    def handle_paper_click(self, item):
        paper = item.data(Qt.UserRole)
        self.update_paper_display(paper)

    def update_paper_display(self, paper):
        self.title_disp.setText(paper.title)
        self.abstract_disp.setText(paper.abstract)

        self.supp_files.clear()
        for f in paper.files:
            list_item = QListWidgetItem()
            custom_item = SuppFileItem(f.url, self, f)
            custom_item.checkbox.setChecked(f.checked)
            list_item.setSizeHint(custom_item.sizeHint())
            self.supp_files.addItem(list_item)
            self.supp_files.setItemWidget(list_item, custom_item)

    def preview_supp_file(self, file_url):
        self.start_load_animation()
        if self.preview_thread.isRunning():
            self.preview_thread.quit()
            self.preview_thread.wait()
        self.preview_thread.file_url = file_url
        self.preview_thread.start()

    def load_preview(self, fname):
        data = self.extract_data_from_file(fname)
        self.update_preview_display(data)

    def extract_data_from_file(self, fname):
        if fname.endswith(('.xlsx', '.xls')):
            data = pd.read_excel(fname, sheet_name=None, header=None, index_col=None)
        elif fname.endswith('.csv'):
            data = pd.read_csv(fname, header=None, index_col=None)
        else: 
            data = pd.read_csv(fname, delimiter='\t', header=None, index_col=None)
        
        os.remove(fname)
        return data

    def update_preview_display(self, data):
        self.stop_load_animation()

        for i in reversed(range(self.previews_layout.count())):
            self.previews_layout.itemAt(i).widget().setParent(None)

        self.determine_data_type_and_display(data)
        self.previews.show()

    def determine_data_type_and_display(self, data):
        if isinstance(data, dict):
            self.display_data_in_tabs(data)
        else:
            self.display_data_in_table(data)

    def display_data_in_tabs(self, data_dict):
        tab_widget = QTabWidget(self.previews)
        for sheet, data in data_dict.items():
            table_widget = QTableWidget()
            self.load_dataframe_to_table(table_widget, data)
            tab_widget.addTab(table_widget, sheet)
        self.previews_layout.addWidget(tab_widget)

    def display_data_in_table(self, data_frame):
        table_widget = QTableWidget(self.previews)
        self.load_dataframe_to_table(table_widget, data_frame)
        self.previews_layout.addWidget(table_widget)

    def load_dataframe_to_table(self, table_widget, data):
        table_widget.setColumnCount(len(data.columns))
        table_widget.setHorizontalHeaderLabels(data.columns.astype(str))
        for row, data in data.iterrows():
            table_widget.insertRow(row)
            for col, value in enumerate(data):
                table_widget.setItem(row, col, QTableWidgetItem(str(value)))

    def on_proceed(self):
        selected_data = self.extract_selected_papers_and_files()
        with open('selected_papers_and_files.json', 'w') as f:
            json.dump(selected_data, f, indent=4)
        QMessageBox.information(self, "Saved", "Selected papers and files have been saved to selected_papers_and_files.json")

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

def main():
    app = QApplication(sys.argv)
    window = CRISPRApp()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()