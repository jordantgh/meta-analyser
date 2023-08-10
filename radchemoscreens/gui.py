import sys
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QListWidget, QListWidgetItem, QTextEdit, QCheckBox, QLabel, QTableWidget, QTableWidgetItem, QProgressBar, QSizePolicy)
from search_for_papers import search_and_summarize_pmc

class SearchThread(QThread):
    article_signal = pyqtSignal(dict, int)  # Will send a single article
    finished_signal = pyqtSignal()  # Signal when all articles are processed

    def __init__(self):
        super().__init__()
        self.query = ""

    def run(self):
        try:
            search_and_summarize_pmc(self.query, callback=self.emit_article)
            self.finished_signal.emit()
        except Exception as e:
            print(e)   # or show a message box, or use any other way to inform the user

    def emit_article(self, article, progress):
        self.article_signal.emit(article, progress)


class CustomListItem(QWidget):
    def __init__(self, title):
        super().__init__()
        self.layout = QHBoxLayout(self)
        self.checkbox = QCheckBox()
        self.label = QLabel(title)
        
        # Set the label's size policy to Fixed
        self.label.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)

        self.layout.addWidget(self.checkbox)
        self.layout.addWidget(self.label)
        self.layout.addStretch(1)  # Add a stretching spacer item
        self.setLayout(self.layout)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)  # Retain this call
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)



class CRISPRApp(QMainWindow):
    def __init__(self):
        super().__init__()

        # Central widget and layout
        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)

        # Left Layout (List of papers)
        self.left_layout = QVBoxLayout()
        self.query_input = QLineEdit(self)
        self.search_status_label = QLabel("")
        self.left_layout.addWidget(self.search_status_label)
        self.search_button = QPushButton("Search", self)
        self.paper_list = QListWidget(self)

        self.left_layout.addWidget(QLabel("Enter Query:"))
        self.left_layout.addWidget(self.query_input)
        self.left_layout.addWidget(self.search_button)
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setRange(0, 100)  # Set range as percentage
        self.progress_bar.setValue(0)
        self.progress_bar.hide()  # Hide initially

        self.left_layout.addWidget(self.progress_bar)
        self.left_layout.addWidget(self.paper_list)

        # Right Layout (Paper details)
        self.right_layout = QVBoxLayout()
        self.title_display = QTextEdit(self)
        self.title_display.setPlaceholderText("Title will be shown here")
        self.abstract_display = QTextEdit(self)
        self.abstract_display.setPlaceholderText("Abstract will be shown here")
        self.supplemental_files_table = QTableWidget(self)
        self.supplemental_files_table.setColumnCount(1)

        self.right_layout.addWidget(QLabel("Title:"))
        self.right_layout.addWidget(self.title_display)
        self.right_layout.addWidget(QLabel("Abstract:"))
        self.right_layout.addWidget(self.abstract_display)
        self.right_layout.addWidget(QLabel("Supplementary Files:"))
        self.right_layout.addWidget(self.supplemental_files_table)

        # Add both layouts to the main layout
        self.main_layout.addLayout(self.left_layout)
        self.main_layout.addLayout(self.right_layout)

        # Connect the signals
        self.search_button.clicked.connect(self.search_papers)
        self.paper_list.itemClicked.connect(self.show_paper_details)

        # Create the search thread and connect its signals
        self.search_thread = SearchThread()
        self.search_thread.article_signal.connect(self.add_article_to_list)
        self.search_thread.finished_signal.connect(self.on_search_finished)
    
    def add_article_to_list(self, article, progress):
        item = QListWidgetItem()
        custom_item = CustomListItem(article["Title"])
        item.setSizeHint(custom_item.sizeHint())
        paper_data = {
            "title": article["Title"],
            "abstract": article["Abstract"],
            "files": []
        }
        item.setData(Qt.UserRole, paper_data)
        self.paper_list.addItem(item)
        self.paper_list.setItemWidget(item, custom_item)
        
        self.progress_bar.setValue(progress)  # Update progress bar value

    def search_papers(self):
        query = self.query_input.text()
        if not query:
            # If the query input is empty, do nothing
            return

        self.paper_list.clear()  # Clear the list before starting a new search
        self.progress_bar.setValue(0)  # Reset progress bar
        self.progress_bar.show()  # Show progress bar
        self.search_status_label.setText("Searching...")  # Start search feedback
        
        # Set the query for the search thread and start the thread
        self.search_thread.query = query
        self.search_thread.start()

    def on_search_finished(self):
        self.progress_bar.hide()  # Hide the progress bar when done
        self.search_status_label.clear()  # Clear search feedback


    def show_paper_details(self, item):
        # Retrieve the paper details from the item's custom data
        paper = item.data(Qt.UserRole)

        self.title_display.setText(paper["title"])
        self.abstract_display.setText(paper["abstract"])

        # Clear and populate supplementary files
        self.supplemental_files_table.setRowCount(0)
        for file in paper["files"]:
            row_position = self.supplemental_files_table.rowCount()
            self.supplemental_files_table.insertRow(row_position)
            self.supplemental_files_table.setItem(row_position, 0, QTableWidgetItem(file))

app = QApplication(sys.argv)
window = CRISPRApp()
window.show()
sys.exit(app.exec_())
