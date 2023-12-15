from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from views.custom_components import TabPage
    from PyQt5.QtGui import QKeyEvent

from PyQt5.QtCore import Qt, QTimer, QObject, pyqtSignal
from PyQt5.QtWidgets import QLabel, QProgressBar, QPushButton, QListWidget, QTabWidget, QTextBrowser, QLineEdit, QSizePolicy


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
        self.outer_tab_widget.addTab(self.metadata_view, "Metadata")

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

    def emit_filter_identity(self):
        self.filter_sig.emit(self.page_identity)

    def emit_prune_identity(self):
        self.prune_sig.emit(self.page_identity)
