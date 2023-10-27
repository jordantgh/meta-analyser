from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from model.article_managers import BaseData, Article, SuppFile, ProcessedTable
    from utils.constants import PageIdentity
    from PyQt5.QtGui import QMouseEvent
    from PyQt5.QtWidgets import QListWidget
    from views.view import View

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLabel, QCheckBox, QSizePolicy
from PyQt5.QtGui import QFontMetrics


class ListItem(QWidget):
    def __init__(
        self,
        data: 'BaseData',
        title: 'str',
        context: 'PageIdentity',
        skip_checkbox=False
    ):
        super().__init__()
        self.data = data
        self.context = context
        self.checkbox = QCheckBox()
        # hack for now, setChecked wants a boolean which Article.checked is not
        if not skip_checkbox:
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
        if self.data.checkbox_togglable:
            self.data.checkbox_toggled()

    def register_observer(self, context):
        self.data.register_observer(self, context)

    def remove(self):
        self.data.remove_observer(self.context)

    def mousePressEvent(self, event: 'QMouseEvent'):
        super().mousePressEvent(event)
        list_widget: 'QListWidget' = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)


class ArticleListItem(ListItem):
    def __init__(self, article: 'Article', context: 'PageIdentity'):
        super().__init__(article, article.title, context, skip_checkbox=True)
        self.register_observer(context)
        self.load_checked_state()

    def load_checked_state(self):
        checked_state = self.data.checked[self.context]
        self.checkbox.setChecked(checked_state)

    def checkbox_toggled(self):
        self.data.checked[self.context] = self.checkbox.isChecked()
        self.data.notify_observers(self.context)

    def update(self, article: 'Article'):
        new_checked_state = article.checked[self.context]
        self.checkbox.toggled.disconnect(self.checkbox_toggled)
        self.checkbox.setChecked(new_checked_state)
        self.checkbox.toggled.connect(self.checkbox_toggled)


class DataListItem(ListItem):
    preview_requested = pyqtSignal(object, object)

    def __init__(
        self,
        main_window: 'View',
        file_data: 'BaseData',
        disp_name_func: 'function',
        context: 'PageIdentity'
    ):
        self.main_window = main_window
        self.page = main_window.active_elements
        self.file_id = file_data.id
        disp_name: 'str' = disp_name_func(file_data)
        super().__init__(file_data, disp_name, context)

    def mousePressEvent(self, event: 'QMouseEvent'):
        super().mousePressEvent(event)
        self.preview_requested.emit(self.data, self.context)

    def get_disp_name(self, text: 'str') -> 'str':
        font_metrics = QFontMetrics(self.main_window.font())
        available_width = self.page.data_list_view.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)


class SuppFileListItem(DataListItem):
    def __init__(
        self,
        main_window: 'View',
        data: 'SuppFile',
        context: 'PageIdentity'
    ):
        super().__init__(
            main_window,
            data,
            lambda fd: self.get_disp_name(fd.url.split('/')[-1]),
            context
        )

        self.file_url = data.url


class ProcessedTableListItem(DataListItem):
    def __init__(
        self,
        main_window: 'View',
        data: 'ProcessedTable',
        context: 'PageIdentity'
    ):
        super().__init__(
            main_window,
            data,
            lambda fd: self.get_disp_name(fd.id),
            context
        )

        self.register_observer(context)

    def checkbox_toggled(self):
        self.data.checked = self.checkbox.isChecked()
        if self.data.checkbox_togglable:
            self.data.checkbox_toggled(self.context)

    def update(self, processed_table: 'ProcessedTable'):
        # Disconnect the signal before updating the checkbox to avoid triggering the signal again, then reconnect it
        self.checkbox.toggled.disconnect(self.checkbox_toggled)
        self.checkbox.setChecked(processed_table.checked)
        self.checkbox.toggled.connect(self.checkbox_toggled)
