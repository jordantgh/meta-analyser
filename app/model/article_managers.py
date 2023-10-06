class BaseData:
    def alert_observers(self):
        return False


class SuppFile(BaseData):
    def __init__(self, article_id, url, id):
        self.checked = True
        self.article_id = article_id
        self.url = url
        self.id = id


class ProcessedTable(BaseData):
    def __init__(self, article_id, id, file_id, num_columns=None):
        self.checked = True
        self.article_id = article_id
        self.id = id
        self.file_id = file_id
        if num_columns is not None:
            self.checked_columns = list(range(num_columns))
        else:
            self.checked_columns = []
        self.observers = []

    def alert_observers(self):
        return True
     
    def register_observer(self, observer):
        self.observers.append(observer)

    def remove_observer(self, observer):
        self.observers.remove(observer)

    def notify_observers(self):
        for observer in self.observers:
            observer.update(self)

    def checkbox_toggled(self):
        self.notify_observers()

    def set_checked(self, state):
        was_checked = self.checked
        self.checked = state
        if state != was_checked:
            self.notify_observers()

class SuppFileManager:
    def __init__(self):
        self.supp_files = {}

    def add_file(self, file):
        self.supp_files[file.id] = file

    def get_file(self, file_id):
        return self.supp_files.get(file_id)
      
    def reset(self):
        self.supp_files = {}


class ProcessedTableManager:
    def __init__(self):
        self.processed_tables = {}

    def add_processed_table(self, table):
        self.processed_tables[table.id] = table

    def get_processed_table(self, table_id):
        return self.processed_tables.get(table_id)
    
    def reset(self):
        self.processed_tables = {}


class Article(BaseData):
    def __init__(
        self, title, abstract, pmc_id, supp_files=[], processed_tables=[]):
        self.checked = True
        self.title = title
        self.abstract = abstract
        self.pmc_id = pmc_id
        self.supp_files = supp_files
        self.processed_tables = processed_tables

    def get_file(self, file_id):
        return next((f for f in self.supp_files if f.id == file_id), None)
      
    def get_table_by_id(self, table_id):
        return next((t for t in self.processed_tables if t.id == table_id), None)


class Bibliography:
    def __init__(self):
        self.articles = {}

    def add_article(self, article):
        self.articles[article.pmc_id] = article

    def get_article(self, article_id):
        return self.articles.get(article_id)

    def get_selected_articles(self):
        return [p for p in self.articles.values() if p.checked]
      
    def reset(self):
        self.articles = {}