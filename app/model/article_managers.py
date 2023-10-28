from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from views.list import ListItem, DataListItem, ArticleListItem

from utils.constants import PageIdentity
from uuid import uuid4, UUID


class BaseData:
    def __init__(self):
        self.checked = True
        self.observers: 'dict[PageIdentity, ListItem]' = {}

    def register_observer(
        self,
        observer: 'ListItem',
        context: 'PageIdentity'
    ):
        self.observers[context] = observer

    def notify_observers(self, context: 'PageIdentity'):
        self.observers[context].update(self)

    def alert_observers(self):
        return True

    def checkbox_togglable(self):
        return True

    def remove_observer(self, context):
        del self.observers[context]

    def stash_observers(self, global_stash):
        if hasattr(self, 'observers'):
            global_stash[id(self)] = self.observers
            self.observers = {}

    def restore_observers(self, global_stash):
        if id(self) in global_stash:
            self.observers = global_stash[id(self)]
            del global_stash[id(self)]


class SuppFile(BaseData):
    def __init__(
            self,
            article: 'Article',
            url: 'str',
            metadata: 'str',
            id: 'UUID'
    ):
        super().__init__()
        self.article = article
        self.article_id = article.pmc_id
        self.url = url
        self.metadata = metadata
        self.id = id

    # don't use the observer pattern for supp files, they don't need it
    def register_observer():
        pass

    def notify_observers():
        pass

    def alert_observers(self):
        return False

    def checkbox_toggled(self):
        self.article.update_based_on_elements(PageIdentity.SEARCH)


class ProcessedTable(BaseData):
    def __init__(
        self, article: 'Article',
        id: 'str',
        file_id: 'UUID',
        num_columns=None
    ):
        super().__init__()
        self.article: 'Article' = article
        self.article_id: 'str' = article.pmc_id
        self.supp_file: 'SuppFile' = article.get_file(file_id)
        self.id: 'str' = id
        self.file_id = file_id
        self.pruned_columns = []
        if num_columns is not None:
            self.checked_columns = list(range(num_columns))
        else:
            self.checked_columns = []
        self.observers: 'dict[PageIdentity, DataListItem]' = {}

    def checkbox_toggled(self, context: 'PageIdentity'):
        self.notify_observers(context)
        self.article.update_based_on_elements(context)

    def set_checked_state(self, checked_state: 'bool', context: 'PageIdentity'):
        old_state = self.checked
        self.checked = checked_state
        if checked_state != old_state:
            self.article.update_based_on_elements(context)
            # first check if observers are registered yet (checking/unchecking
            # can happen before the table is displayed via the filter method)
            if context in self.observers:
                self.notify_observers(context)


class ProcessedTableManager:
    def __init__(self):
        self.processed_tables = {}

    def add_processed_table(self, table: 'ProcessedTable'):
        self.processed_tables[table.id] = table

    def get_processed_table(self, table_id: 'str') -> 'ProcessedTable':
        return self.processed_tables.get(table_id)

    def reset(self):
        self.processed_tables = {}


class Article(BaseData):
    def __init__(self, article_json: 'dict'):
        super().__init__()
        self.checked: 'dict[PageIdentity, bool]' = {
            context: True for context in PageIdentity
        }

        self.title: 'str' = article_json["Title"]
        self.authors: 'str' = article_json["Authors"]
        self.abstract: 'str' = article_json["Abstract"]
        self.pmc_id: 'str' = article_json["PMCID"]
        self.url: 'str' = article_json["URL"]

        supp_data: 'dict[str, str]' = article_json["SupplementaryFiles"]

        self.supp_files = [
            SuppFile(self, url, meta, uuid4())
            for url, meta in supp_data.items()
        ]

        self.processed_tables: 'list[ProcessedTable]' = []
        self.pruned_tables: 'list[ProcessedTable]' = []
        self.observers: 'dict[PageIdentity, ArticleListItem]' = {}

    def cascade_checked_state(
            self,
            context: 'PageIdentity',
            is_checked: 'bool' = None
    ):
        if is_checked is None:
            is_checked = self.checked[context]

        hierarchy = [context for context in PageIdentity]
        current_index = hierarchy.index(context)
        self.checked[context] = is_checked

        if current_index == len(hierarchy) - 1:
            return

        for sub_context in hierarchy[current_index + 1:]:
            self.cascade_checked_state(sub_context, is_checked)

    def update_based_on_elements(self, context: 'PageIdentity'):
        has_checked = self.has_checked_elements(context)
        self.checked[context] = has_checked
        self.notify_observers(context)

    def has_checked_elements(self, context: 'PageIdentity') -> 'bool':
        if context == PageIdentity.SEARCH:
            return any(f.checked for f in self.supp_files)
        elif context == PageIdentity.PARSED:
            return any(t.checked for t in self.processed_tables)
        elif context == PageIdentity.PRUNED:
            return any(t.checked for t in self.pruned_tables)

    def get_file(self, file_id: 'UUID') -> 'SuppFile':
        return next((f for f in self.supp_files if f.id == file_id), None)

    def get_table_by_id(self, table_id: 'str') -> 'ProcessedTable':
        return next(
            (t for t in self.processed_tables if t.id == table_id), None
        )


class Bibliography:
    def __init__(self):
        self.articles = {}

    def add_article(self, article: 'Article'):
        self.articles[article.pmc_id] = article

    def get_article(self, article_id: 'str') -> 'Article':
        return self.articles.get(article_id)

    def get_selected_articles(self, context: 'PageIdentity') -> 'list[Article]':
        selected = []
        article: 'Article'
        for article in self.articles.values():
            if article.checked[context]:
                selected.append(article)

        return selected

    def reset(self):
        self.articles = {}


def stash_all_observers(root_object, global_stash, visited_objects):
    object_id = id(root_object)
    if object_id in visited_objects:
        return
    visited_objects.add(object_id)

    if hasattr(root_object, 'stash_observers'):
        root_object.stash_observers(global_stash)

    for attr_name in dir(root_object):
        attr = getattr(root_object, attr_name)
        if isinstance(attr, list):
            for item in attr:
                stash_all_observers(item, global_stash, visited_objects)
        elif isinstance(attr, dict):
            for item in attr.values():
                stash_all_observers(item, global_stash, visited_objects)
        elif isinstance(attr, BaseData):
            stash_all_observers(attr, global_stash, visited_objects)


def restore_all_observers(root_object, global_stash, visited_objects):
    object_id = id(root_object)
    if object_id in visited_objects:
        return
    visited_objects.add(object_id)

    if hasattr(root_object, 'restore_observers'):
        root_object.restore_observers(global_stash)

    for attr_name in dir(root_object):
        attr = getattr(root_object, attr_name)
        if isinstance(attr, list):
            for item in attr:
                restore_all_observers(item, global_stash, visited_objects)
        elif isinstance(attr, dict):
            for item in attr.values():
                restore_all_observers(item, global_stash, visited_objects)
        elif isinstance(attr, BaseData):
            restore_all_observers(attr, global_stash, visited_objects)
