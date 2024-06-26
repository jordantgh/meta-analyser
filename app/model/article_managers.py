from typing import TYPE_CHECKING, Optional
if TYPE_CHECKING:
    from views.list import ListItem, DataListItem, ArticleListItem
    from model.database import TableDBManager

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


class ColMapping:
    def __init__(self, table: 'ProcessedTable', colour: 'Optional[str]' = None):
        self.colour = colour
        self.table = table
        self.ascending = True
        self.ids: 'list[tuple]' = []
        self.values: 'list[tuple]' = []

    def add_id_col(self, id_col: 'tuple'):
        self.ids.append(id_col)

    def add_value_col(self, value_col: 'tuple'):
        self.values.append(value_col)

    def set_order(self, ascending: 'bool'):
        self.ascending = ascending


class ProcessedTable(BaseData):
    def __init__(
        self,
        article: 'Article',
        id: 'str',
        file_id: 'UUID',
        num_columns: 'int' = None
    ):
        super().__init__()
        self.article = article
        self.id = id
        self.file_id = file_id
        self.highlighting_enabled = False
        self.mappings: 'list[ColMapping]' = []

        if num_columns is not None:
            self.checked_columns: 'list[int]' = list(range(num_columns))
        else:
            self.checked_columns: 'list[int]' = []

        self.article_id: 'str' = article.pmc_id
        self.supp_file: 'SuppFile' = article.get_file(file_id)
        self.pruned_columns: 'list[int]' = []
        self.observers: 'dict[PageIdentity, DataListItem]' = {}
        self.tags = []

    def add_mapping(self, colour: 'Optional[str]' = None):
        mapping = ColMapping(self, colour)
        self.mappings.append(mapping)

    def get_existing_mapping_colours(self):
        return [mapping.colour for mapping in self.mappings]

    def clear_mappings(self):
        self.mappings.clear()

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

    def add_tag(self, tag: 'str', db_manager: 'TableDBManager'):
        if tag not in self.tags:
            self.tags.append(tag)
            db_manager.update_table_tags(self.id, self.tags)

    def remove_tag(self, tag: 'str', db_manager: 'TableDBManager'):
        if tag in self.tags:
            self.tags.remove(tag)
            db_manager.update_table_tags(self.id, self.tags)

    def get_tags(self):
        return self.tags


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

        self.checked: 'dict[PageIdentity, bool]' = {
            context: True for context in PageIdentity
        }

        self.processed_tables: 'list[ProcessedTable]' = []
        self.pruned_tables: 'list[ProcessedTable]' = []
        self.observers: 'dict[PageIdentity, ArticleListItem]' = {}

    def checkbox_toggled(self, context: 'PageIdentity'):
        # ensure this only works if the checked state would be set to False
        # (to prevent all tables being selected when the article is selected)
        if self.checked[context] != True:
            # cascade state change down to the tables
            new_checked_state = self.checked[context]
            tables = self._get_tables_by_context(context)
            for table in tables:
                table.set_checked_state(new_checked_state, context)

    def _get_tables_by_context(
        self, context: 'PageIdentity'
    ) -> 'list[ProcessedTable]':
        # Helper method to get tables by context
        if context == PageIdentity.PARSED:
            return self.processed_tables
        if context == PageIdentity.PRUNED:
            return self.pruned_tables
        return []

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
        for con in (PageIdentity.PARSED, PageIdentity.PRUNED):
            # check the context exists in the observers dict first
            # this is mainly for when there hasnt yet been any pruning

            # Note that iterating through article contexts like this is a
            # hack until I figure out whether I need page specific checked
            # states at all. (TODO)

            if con in self.observers:
                self.checked[con] = has_checked
                self.notify_observers(con)

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


# Helpers for saving
# UI elements are not serializable, so need to be stashed before saving

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
