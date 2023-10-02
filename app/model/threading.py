from PyQt5.QtCore import QThread, pyqtSignal

from model.file_io import download_supp, extract_dfs
from model.tabular_operations import parse_tables
from scripts.search_for_papers import query_pmc

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
        query_pmc(self.query, callback=self.article_sig.emit, thread = self)
        self.finished_sig.emit()


class FilePreviewThread(QThread):
    prev_ready_sig = pyqtSignal(dict)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url

    def run(self):
        fname = download_supp(self.file_url)
        data = extract_dfs(fname)
        self.prev_ready_sig.emit(data)


class FileProcessingThread(QThread):
    article_sig = pyqtSignal(object, list, int)
    finished_sig = pyqtSignal()

    def __init__(self, db_manager):
        super().__init__()
        self.selected_articles = []
        self.db_manager = db_manager

    def run(self):
        parse_tables(self.selected_articles, self.db_manager, callback=self.article_sig.emit)
        self.finished_sig.emit()