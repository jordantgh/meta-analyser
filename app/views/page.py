from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QLabel, QProgressBar, QPushButton, QListWidget, QTextEdit, QLineEdit, QWidget, QVBoxLayout

class CommonPageElements:
    def __init__(self, page):
        self.prog_bar = QProgressBar(page)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()
        self.article_list = QListWidget(page)
        self.title_disp = QTextEdit(page)
        self.title_disp.setPlaceholderText("Title will be shown here")
        self.abstract_disp = QTextEdit(page)
        self.abstract_disp.setPlaceholderText("Abstract will be shown here")
        self.supp_files_view = QListWidget(page)
        self.previews = QWidget(page)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)
        self.loading_label = QLabel(page)
        self.loading_label.setAlignment(Qt.AlignCenter)

class SearchPageElements(CommonPageElements):
    def __init__(self, page):
        super().__init__(page)
        self.search_status = QLabel(page)
        self.query_field = QLineEdit(page)
        self.search_btn = QPushButton("Search", page)
        self.stop_search_btn = QPushButton("Stop Search", page)
        self.stop_search_btn.hide()
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", page)
        

class ProcessedPageElements(CommonPageElements):
    def __init__(self, page):
        super().__init__(page)
        self.query_filter_field = QLineEdit(page)
        self.filter_btn = QPushButton("Filter", page)
        self.prune_btn = QPushButton("Prune Tables and Columns", page)