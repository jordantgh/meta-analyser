import sys
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import *
from search_for_papers import query_pmc
import logging

class SearchThread(QThread):
    article_signal = pyqtSignal(dict, int)  # Will send a single article
    finished_signal = pyqtSignal()  # Signal when all articles are processed

    def __init__(self):
        super().__init__()
        self.query = ""

    def run(self):
        try:
            query_pmc(self.query, callback=self.emit_article)
            self.finished_signal.emit()
        except Exception as e:
            logging.error(f"Unhandled exception in SearchThread: {e}")
            print(e)

    def emit_article(self, article, progress):
        self.article_signal.emit(article, progress)


class UIListItem(QWidget):
    def __init__(self, title):
        super().__init__()
        layout = QHBoxLayout(self)
        checkbox = QCheckBox()
        label = QLabel(title)
        
        # Set the label's size policy to Fixed
        label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        layout.addWidget(checkbox)
        layout.addWidget(label)
        layout.addStretch(1)  # Add a stretching spacer to left-align
        self.setLayout(layout)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)

class CRISPRApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Central widget and layout
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        # Left Layout
        left_panel = QVBoxLayout()
        self.query_input = QLineEdit(self)
        self.search_status_label = QLabel("")
        left_panel.addWidget(self.search_status_label)
        self.search_button = QPushButton("Search", self)
        self.paper_list = QListWidget(self)
        left_panel.addWidget(QLabel("Enter Query:"))
        left_panel.addWidget(self.query_input)
        left_panel.addWidget(self.search_button)
        self.prog_bar = QProgressBar(self)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()
        left_panel.addWidget(self.prog_bar)
        left_panel.addWidget(self.paper_list)

        # Right Layout
        right_panel = QVBoxLayout()
        self.title_display = QTextEdit(self)
        self.title_display.setPlaceholderText("Title will be shown here")
        self.abstract_display = QTextEdit(self)
        self.abstract_display.setPlaceholderText("Abstract will be shown here")
        self.supp_files_table = QTableWidget(self)
        self.supp_files_table.setColumnCount(1)
        right_panel.addWidget(QLabel("Title:"))
        right_panel.addWidget(self.title_display)
        right_panel.addWidget(QLabel("Abstract:"))
        right_panel.addWidget(self.abstract_display)
        right_panel.addWidget(QLabel("Supplementary Files:"))
        right_panel.addWidget(self.supp_files_table)

        main_panel = QHBoxLayout(central_widget)
        main_panel.addLayout(left_panel)
        main_panel.addLayout(right_panel)

        # Connect the signals
        self.search_button.clicked.connect(self.search_papers)
        self.paper_list.itemClicked.connect(self.show_paper_details)

        # Create the search thread and connect its signals
        self.search_thread = SearchThread()
        self.search_thread.article_signal.connect(self.add_article_to_list)
        self.search_thread.finished_signal.connect(self.on_search_finished)
    
    def add_article_to_list(self, article, progress):
        item = QListWidgetItem()
        custom_item = UIListItem(article["Title"])
        item.setSizeHint(custom_item.sizeHint())
        paper_data = {
            "title": article["Title"],
            "abstract": article["Abstract"],
            "files": article["SupplementaryFiles"]
        }
        item.setData(Qt.UserRole, paper_data)
        self.paper_list.addItem(item)
        self.paper_list.setItemWidget(item, custom_item)
        self.prog_bar.setValue(progress)

    def search_papers(self):
        query = self.query_input.text()
        if not query: return
        self.paper_list.clear()
        self.prog_bar.setValue(0)
        self.prog_bar.show()
        self.search_status_label.setText("Searching...")
        self.search_thread.query = query
        self.search_thread.start()

    def on_search_finished(self):
        self.prog_bar.hide()
        self.search_status_label.clear()

    def show_paper_details(self, item):
        # Retrieve the paper details from the item's custom data
        paper = item.data(Qt.UserRole)

        self.title_display.setText(paper["title"])
        self.abstract_display.setText(paper["abstract"])

        # Clear and populate supplementary files
        self.supp_files_table.setRowCount(0)
        for file_url in paper["files"]:
            filename = file_url.split('/')[-1]
            row_position = self.supp_files_table.rowCount()
            self.supp_files_table.insertRow(row_position)
            self.supp_files_table.setItem(row_position, 0, QTableWidgetItem(filename))

app = QApplication(sys.argv)
window = CRISPRApp()
window.show()
sys.exit(app.exec_())
