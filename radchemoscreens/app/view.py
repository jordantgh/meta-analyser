from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QPushButton, QListWidget, QListWidgetItem, QProgressBar, QTextEdit, QCheckBox, QPushButton, QTableWidget, QTableWidgetItem, QSizePolicy, QTabWidget

class UIListItem(QWidget):
    def __init__(self, title, data=None):
        super().__init__()
        self.data = data
        self.checkbox = QCheckBox()
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

    # override default mousePressEvent behaviour
    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)

class SuppFileItem(UIListItem):
    preview_requested = pyqtSignal(str)
    def __init__(self, file_url, main_window, file_instance):
        self.file_url = file_url
        self.main_window = main_window
        disp_name = self.get_disp_name(file_url.split('/')[-1])

        super().__init__(disp_name, file_instance)

        self.preview_btn = QPushButton("Preview", self)
        self.preview_btn.clicked.connect(self.preview_file)
        self.layout().addWidget(self.preview_btn)

    def get_disp_name(self, text):
        font_metrics = QFontMetrics(self.main_window.font())
        available_width = self.main_window.supp_files.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)

    def preview_file(self):
        self.preview_requested.emit(self.file_url)

class CRISPRView(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        self.init_uis()
        self.setup_layouts()

    def init_uis(self):
        # Central widget
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        self.stop_search_btn = QPushButton("Stop Search", self)
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", self)

        # Paper list elements
        self.query_field = QLineEdit(self)
        self.search_status = QLabel("")
        self.search_btn = QPushButton("Search", self)
        self.paper_list = QListWidget(self)
        self.prog_bar = QProgressBar(self)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()

        # Paper details elements
        self.title_disp = QTextEdit(self)
        self.title_disp.setPlaceholderText("Title will be shown here")
        self.abstract_disp = QTextEdit(self)
        self.abstract_disp.setPlaceholderText("Abstract will be shown here")
        self.supp_files = QListWidget(self)
        self.load_label = QLabel(self)
        self.load_label.setAlignment(Qt.AlignCenter)
        self.previews = QWidget(self)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)

    def setup_layouts(self):
        # Paper list layout
        pane_0 = QVBoxLayout()
        pane_0.addWidget(QLabel("Enter Query:"))
        pane_0.addWidget(self.query_field)
        pane_0.addWidget(self.prog_bar)
        pane_0.addWidget(self.search_btn)
        pane_0.addWidget(self.stop_search_btn)
        pane_0.addWidget(self.paper_list)

        # Paper details layout
        pane_1 = QVBoxLayout()
        pane_1.addWidget(QLabel("Title:"))
        pane_1.addWidget(self.title_disp)
        pane_1.addWidget(QLabel("Abstract:"))
        pane_1.addWidget(self.abstract_disp)
        pane_1.addWidget(QLabel("Supplementary Files:"))
        pane_1.addWidget(self.supp_files)
        pane_1.addWidget(self.load_label)

        # Main layout
        main_pane = QHBoxLayout(self.centralWidget())
        main_pane.addLayout(pane_0)
        main_pane.addLayout(pane_1)
        main_pane.addWidget(self.previews)
        
    def init_load_animation(self):
        self.load_timer = QTimer(self)
        self.load_dots = 0
        self.load_timer.timeout.connect(self.update_load_text)

    def start_load_animation(self):
        self.load_dots = 0
        self.load_label.setText("Loading.")
        self.load_timer.start(500)  # ms

    def stop_load_animation(self):
        self.load_timer.stop()
        self.load_label.clear()

    def update_load_text(self):
        self.load_dots = (self.load_dots + 1) % 4
        self.load_label.setText("Loading" + "." * self.load_dots)

    def display_paper_in_ui(self, paper_data, progress):
        item = QListWidgetItem()
        paper_name = UIListItem(paper_data.title, paper_data)
        item.setSizeHint(paper_name.sizeHint())
        item.setData(Qt.UserRole, paper_data)
        self.paper_list.addItem(item)
        self.paper_list.setItemWidget(item, paper_name)
        self.prog_bar.setValue(progress + 1)

    def update_paper_display(self, paper):
        self.title_disp.setText(paper.title)
        self.abstract_disp.setText(paper.abstract)

        self.supp_files.clear()
        for f in paper.files:
            list_item = QListWidgetItem()
            custom_item = SuppFileItem(f.url, self, f)
            custom_item.checkbox.setChecked(f.checked)
            list_item.setSizeHint(custom_item.sizeHint())
            self.supp_files.addItem(list_item)
            self.supp_files.setItemWidget(list_item, custom_item)

    def display_data_in_tabs(self, data_dict):
        tab_widget = QTabWidget(self.previews)
        for sheet, data in data_dict.items():
            table_widget = QTableWidget()
            self.load_dataframe_to_table(table_widget, data)
            tab_widget.addTab(table_widget, sheet)
        self.previews_layout.addWidget(tab_widget)

    def display_data_in_table(self, data_frame):
        table_widget = QTableWidget(self.previews)
        self.load_dataframe_to_table(table_widget, data_frame)
        self.previews_layout.addWidget(table_widget)

    def load_dataframe_to_table(self, table_widget, data):
        table_widget.setColumnCount(len(data.columns))
        table_widget.setHorizontalHeaderLabels(data.columns.astype(str))
        for row, data in data.iterrows():
            table_widget.insertRow(row)
            for col, value in enumerate(data):
                table_widget.setItem(row, col, QTableWidgetItem(str(value)))