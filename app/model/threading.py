from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from model.database import TableDBManager
    from model.article_managers import Article

from PyQt5.QtCore import QThread, pyqtSignal

from model.file_io import download_supp, extract_dfs
from model.tabular_operations import parse_tables
from scripts.search_for_papers import query_pmc


class SearchThread(QThread):
    article_sig = pyqtSignal(dict, int)
    finished_sig = pyqtSignal(object)

    def __init__(self):
        super().__init__()
        self.query = ""
        self.should_stop = False

    def stop(self):
        self.should_stop = True

    def prepare(self, query: 'str'):
        self.query = query
        self.should_stop = False

    def run(self):
        query_pmc(self.query, callback=self.article_sig.emit, thread=self)
        self.finished_sig.emit(self)


class FilePreviewThread(QThread):
    prev_ready_sig = pyqtSignal(dict)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url
        self.should_stop = False

    def stop(self):
        self.should_stop = True

    def prepare(self, file_url: 'str'):
        self.file_url = file_url
        self.should_stop = False

    def run(self):
        fname = download_supp(self.file_url, self.should_stop)
        if fname is None:
            return
        data = extract_dfs(fname, self.should_stop)
        if data is not None:
            self.prev_ready_sig.emit(data)


class FileProcessingThread(QThread):
    article_sig = pyqtSignal(object, int, list)
    finished_sig = pyqtSignal()

    def __init__(self, db_manager: 'TableDBManager'):
        super().__init__()
        self.selected_articles = []
        self.db_manager = db_manager
        self.should_stop = False

    def prepare(self, selected_articles: 'list[Article]'):
        self.selected_articles = selected_articles
        self.should_stop = False

    def run(self):
        parse_tables(
            self.selected_articles,
            self.db_manager,
            self.should_stop,
            callback=self.article_sig.emit
        )

        self.finished_sig.emit()
