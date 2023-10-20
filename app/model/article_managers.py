from utils.constants import PageIdentity


class BaseData:
    def alert_observers(self):
        return False

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
    def __init__(self, article, url, id):
        self.checked = True
        self.article = article
        self.article_id = article.pmc_id
        self.url = url
        self.id = id

    def checkbox_toggled(self):
        self.article.update_based_on_elements(PageIdentity.SEARCH)


class ProcessedTable(BaseData):
    def __init__(self, article, id, file_id, num_columns=None):
        self.checked = True
        self.article = article
        self.article_id = article.pmc_id
        self.id = id
        self.file_id = file_id
        self.pruned_columns = []
        if num_columns is not None:
            self.checked_columns = list(range(num_columns))
        else:
            self.checked_columns = []
        self.observers = {}

    def register_observer(self, observer, context):
        self.observers[context] = observer

    def notify_observers(self, context):
        self.observers[context].update(self)

    def alert_observers(self):
        return True

    def checkbox_toggled(self, context):
        self.notify_observers(context)
        self.article.update_based_on_elements(context)

    def set_checked(self, state, context):
        was_checked = self.checked
        self.checked = state
        if state != was_checked:
            self.article.update_based_on_elements(context)
            # first check if observers are registered yet (checking/unchecking
            # can happen before the table is displayed via the filter method)
            if context in self.observers:
                self.notify_observers(context)


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
            self,
            title,
            authors,
            abstract,
            pmc_id,
            url,
            supp_files=[],
            processed_tables=[]
    ):
        self.checked = {context: True for context in PageIdentity}
        self.title = title
        self.authors = authors
        self.abstract = abstract
        self.pmc_id = pmc_id
        self.url = url
        self.supp_files = supp_files
        self.processed_tables = processed_tables
        self.pruned_tables = []
        self.observers = {}

    def cascade_checked_state(self, context, is_checked=None):
        if is_checked is None:
            is_checked = self.checked[context]

        hierarchy = [context for context in PageIdentity]

        # Find the position of the current context in the hierarchy
        index = hierarchy.index(context)

        # Update the checked state of the parent context (preserves initial
        # state for the initially called context)
        self.checked[context] = is_checked

        # Stop recursion if we're at the last element in the hierarchy
        if index == len(hierarchy) - 1:
            return

        for sub_context in hierarchy[index + 1:]:
            self.cascade_checked_state(sub_context, is_checked)

    def alert_observers(self):
        return True

    def register_observer(self, observer, context):
        self.observers[context] = observer

    def notify_observers(self, context):
        self.observers[context].update(self)

    def update_based_on_elements(self, context):
        has_checked = self.has_checked_elements(context)
        self.checked[context] = has_checked
        self.notify_observers(context)

    def has_checked_elements(self, context):
        if context == PageIdentity.SEARCH:
            return any(f.checked for f in self.supp_files)
        elif context == PageIdentity.PARSED:
            return any(t.checked for t in self.processed_tables)
        elif context == PageIdentity.PRUNED:
            return any(t.checked for t in self.pruned_tables)

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

    def get_selected_articles(self, context):
        selected = []
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
