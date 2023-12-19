from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from views.custom_components import TabPage
    from PyQt5.QtGui import QKeyEvent

from PyQt5.QtCore import Qt, QTimer, QObject, pyqtSignal
from PyQt5.QtWidgets import QLabel, QProgressBar, QPushButton, QListWidget, QTabWidget, QTextBrowser, QLineEdit, QSizePolicy, QHBoxLayout, QVBoxLayout, QWidget, QListWidgetItem


class QPushButton(QPushButton):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def keyPressEvent(self, event: 'QKeyEvent'):
        if event is not None:
            if event.key() == Qt.Key_Enter or event.key() == Qt.Key_Return:
                self.setDown(True)
                QTimer.singleShot(100, self.simulateClick)
                return
        super().keyPressEvent(event)

    def simulateClick(self):
        self.click()
        self.setDown(False)


class QListWidget(QListWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def keyPressEvent(self, event: 'QKeyEvent'):
        super().keyPressEvent(event)
        if event.key() in [Qt.Key_Up, Qt.Key_Down]:
            self.itemClicked.emit(self.currentItem())


class TagEntryWidget(QLineEdit):
    tagAdded = pyqtSignal(str)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setPlaceholderText("Type a descriptive tag and press Enter...")

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            text = self.text().strip().replace(" ", "_")
            if text:
                self.tagAdded.emit(text)
                self.clear()
        super().keyPressEvent(event)


class TagWidget(QWidget):
    # Custom widget that includes a label and a delete button
    removeTag = pyqtSignal(str)

    def __init__(self, tag, parent=None):
        super().__init__(parent)
        self.tag = tag

        # Layout
        layout = QHBoxLayout(self)
        tag_label = QLabel(f"#{tag}")
        remove_button = QPushButton("X")
        remove_button.setFixedSize(20, 20)  # Small, fixed-size button
        layout.addWidget(tag_label)
        layout.addWidget(remove_button)
        layout.addStretch()
        layout.setContentsMargins(0, 0, 0, 0)

        # Connect the remove button signal
        remove_button.clicked.connect(self.emitRemoveSignal)

    def emitRemoveSignal(self):
        self.removeTag.emit(self.tag)

class TagsDisplayWidget(QListWidget):
    tagRemoved = pyqtSignal(str)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setSelectionMode(QListWidget.NoSelection)

    def addTag(self, tag):
        item = QListWidgetItem(self)
        tag_widget = TagWidget(tag, self)
        item.setSizeHint(tag_widget.sizeHint())
        self.addItem(item)
        self.setItemWidget(item, tag_widget)
        tag_widget.removeTag.connect(self.removeTag)

    def removeTag(self, tag):
        # Find the item and remove it from the list
        for index in range(self.count()):
            item = self.item(index)
            widget = self.itemWidget(item)
            if widget.tag == tag:
                self.takeItem(index)
                break
        self.tagRemoved.emit(tag)

class PageElements(QObject):  # QObject needed for signalling
    def __init__(self, parent_tab: 'TabPage'):
        super().__init__()
        self.page_identity = parent_tab.page_identity
        self.prog_bar = QProgressBar(parent_tab)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()
        self.article_ui_list = QListWidget(parent_tab)
        self.data_ui_list = QListWidget(parent_tab)

        self.title_abstract_disp = QTextBrowser(parent_tab)
        self.title_abstract_disp.setMinimumHeight(100)
        self.title_abstract_disp.setPlaceholderText(
            "Title/abstract will be shown here"
        )

        self.title_abstract_disp.setOpenExternalLinks(True)
        self.title_abstract_disp.setSizePolicy(
            QSizePolicy.Expanding, QSizePolicy.Preferred
        )

        self.title_abstract_disp.setFocusPolicy(Qt.NoFocus)

        self.data_previews = QTabWidget()

        self.metadata_view = QTextBrowser()
        self.metadata_view.setOpenExternalLinks(True)

        self.info_tabs_widget = QTabWidget(parent_tab)
        self.info_tabs_widget.addTab(self.data_previews, "Previews")
        self.info_tabs_widget.addTab(self.metadata_view, "Metadata")

        self.loading_label = QLabel(parent_tab)
        self.loading_label.setAlignment(Qt.AlignCenter)


class SearchPageElements(PageElements):
    def __init__(self, parent_tab: 'TabPage'):
        super().__init__(parent_tab)
        self.search_status = QLabel(parent_tab)
        self.query_field = QLineEdit(parent_tab)
        self.search_btn = QPushButton("Search", parent_tab)
        self.query_field.returnPressed.connect(self.search_btn.click)
        self.stop_search_btn = QPushButton("Stop Search", parent_tab)
        self.stop_search_btn.hide()
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", parent_tab)


class ProcessedPageElements(PageElements):
    filter_sig = pyqtSignal(object)
    prune_sig = pyqtSignal(object)

    def __init__(self, parent_tab: 'TabPage'):
        super().__init__(parent_tab)
        self.query_filter_field = QLineEdit(parent_tab)
        self.filter_btn = QPushButton("Filter", parent_tab)
        self.query_filter_field.returnPressed.connect(self.filter_btn.click)
        self.filter_btn.clicked.connect(self.emit_filter_identity)

        self.prune_btn = QPushButton("Prune Tables and Columns", parent_tab)
        self.prune_btn.clicked.connect(self.emit_prune_identity)

        self.tags_entry_widget = TagEntryWidget(parent_tab)
        self.tags_display_widget = TagsDisplayWidget(parent_tab)
        self.tags_layout = QHBoxLayout()
        self.tags_layout.addWidget(self.tags_entry_widget)
        self.tags_display_layout = QVBoxLayout()
        self.tags_display_layout.addWidget(self.tags_display_widget)

        self.tags_widget = QWidget()
        self.tags_full_layout = QVBoxLayout(self.tags_widget)
        self.tags_full_layout.addLayout(self.tags_layout)
        self.tags_full_layout.addLayout(self.tags_display_layout)

        self.info_tabs_widget.addTab(self.tags_widget, "Tags")

    def emit_filter_identity(self):
        self.filter_sig.emit(self.page_identity)

    def emit_prune_identity(self):
        self.prune_sig.emit(self.page_identity)
