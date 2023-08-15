import sys
from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import *
from search_for_papers import query_pmc
import logging
import requests
import os
import pandas as pd

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

class FilePreviewThread(QThread):
    preview_ready_signal = pyqtSignal(str)

    def __init__(self, file_url):
        super().__init__()
        self.file_url = file_url

    def run(self):
        filename = download_file(self.file_url)
        self.preview_ready_signal.emit(filename)

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

class SuppFileItem(UIListItem):
    def __init__(self, file_url, main_window):
        self.file_url = file_url
        self.main_window = main_window
        
        # Get the width of the parent container (supp_files_list)
        parent_width = self.main_window.supp_files_list.width()

        # Space for preview button, checkbox, and layout margins/paddings
        reserved_width = 150
        
        available_width = parent_width - reserved_width

        filename = file_url.split('/')[-1]
        font_metrics = QFontMetrics(self.main_window.font())

        # Elide the text if it's too long for the available space
        elided_filename = font_metrics.elidedText(filename, Qt.ElideMiddle, available_width)

        super().__init__(elided_filename)

        # Set the size policy of the label
        self.layout().itemAt(1).widget().setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        self.layout().itemAt(1).widget().setMaximumWidth(available_width)
        
        self.preview_btn = QPushButton("Preview", self)
        self.preview_btn.clicked.connect(self.preview_file)
        self.layout().addWidget(self.preview_btn)
    
    def preview_file(self):
        # Use the main_window reference to access the preview_supp_file method
        self.main_window.preview_supp_file(self.file_url)

class CRISPRApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)

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
        self.supp_files_list = QListWidget(self)
        right_panel.addWidget(QLabel("Title:"))
        right_panel.addWidget(self.title_display)
        right_panel.addWidget(QLabel("Abstract:"))
        right_panel.addWidget(self.abstract_display)
        right_panel.addWidget(QLabel("Supplementary Files:"))
        right_panel.addWidget(self.supp_files_list)

        # Add the preview pane to the right
        self.preview_pane = QWidget(self)
        self.preview_pane_layout = QVBoxLayout(self.preview_pane)
        self.preview_pane.setLayout(self.preview_pane_layout)

        main_panel = QHBoxLayout(central_widget)
        main_panel.addLayout(left_panel)
        main_panel.addLayout(right_panel)
        main_panel.addWidget(self.preview_pane)
        self.preview_pane.hide()

        # Connect the signals
        self.search_button.clicked.connect(self.search_papers)
        self.paper_list.itemClicked.connect(self.show_paper_details)

        # Create the search thread and connect its signals
        self.search_thread = SearchThread()
        self.search_thread.article_signal.connect(self.add_article_to_list)
        self.search_thread.finished_signal.connect(self.on_search_finished)

        self.file_preview_thread = FilePreviewThread("")
        self.file_preview_thread.preview_ready_signal.connect(self.load_preview)
    
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
        self.prog_bar.setValue(progress + 1)

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
        self.supp_files_list.clear()  # This replaces the supp files table
        for file_url in paper["files"]:
            list_item = QListWidgetItem()
            custom_item = SuppFileItem(file_url, self)
            list_item.setSizeHint(custom_item.sizeHint())
            self.supp_files_list.addItem(list_item)
            self.supp_files_list.setItemWidget(list_item, custom_item)

    def preview_supp_file(self, file_url):
        # Stop the thread if it's running
        if self.file_preview_thread.isRunning():
            self.file_preview_thread.terminate()
            self.file_preview_thread.wait()

        # Set the file URL and start the thread
        self.file_preview_thread.file_url = file_url
        self.file_preview_thread.start()

    # This method goes inside the CRISPRApp class
    def load_preview(self, filename):
        # Clear any previous content
        for i in reversed(range(self.preview_pane_layout.count())):
            self.preview_pane_layout.itemAt(i).widget().setParent(None)
        
        # Load and display the file content based on its type
        if (filename.endswith('.xlsx') or filename.endswith('.xls')):
            # For Excel files, load all sheets
            data = pd.read_excel(filename, sheet_name=None)
            if isinstance(data, dict):
                # Multiple sheets
                tab_widget = QTabWidget(self.preview_pane)
                for sheet_name, sheet_data in data.items():
                    table_widget = QTableWidget()
                    self.load_dataframe_to_table(table_widget, sheet_data)
                    tab_widget.addTab(table_widget, sheet_name)
                self.preview_pane_layout.addWidget(tab_widget)
            else:
                # Single sheet
                table_widget = QTableWidget(self.preview_pane)
                self.load_dataframe_to_table(table_widget, data)
                self.preview_pane_layout.addWidget(table_widget)
        elif filename.endswith('.csv'):
            data = pd.read_csv(filename)
            table_widget = QTableWidget(self.preview_pane)
            self.load_dataframe_to_table(table_widget, data)
            self.preview_pane_layout.addWidget(table_widget)
        elif filename.endswith('.txt'):
            data = pd.read_csv(filename, delimiter='\t')
            table_widget = QTableWidget(self.preview_pane)
            self.load_dataframe_to_table(table_widget, data)
            self.preview_pane_layout.addWidget(table_widget)

        self.preview_pane.show()
        # Remove the downloaded file after previewing
        os.remove(filename)

    def load_dataframe_to_table(self, table_widget, data):
        # Display the DataFrame in QTableWidget
        table_widget.setColumnCount(len(data.columns))
        table_widget.setHorizontalHeaderLabels(data.columns)
        for row_num, row_data in data.iterrows():
            table_widget.insertRow(row_num)
            for col_num, value in enumerate(row_data):
                table_widget.setItem(row_num, col_num, QTableWidgetItem(str(value)))

def download_file(url):
    local_filename = url.split('/')[-1]
    if not local_filename.endswith(('.txt', '.csv', '.xlsx')):
        print(f"Skipping {local_filename} because it's not a .txt, .csv or .xlsx file.")
        return
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

app = QApplication(sys.argv)
window = CRISPRApp()
window.show()
sys.exit(app.exec_())
