from PyQt5.QtGui import QFontMetrics, QStandardItemModel, QStandardItem
from PyQt5.QtCore import Qt, pyqtSignal, QTimer
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QPushButton, QListWidget, \
    QListWidgetItem, QProgressBar, QTextEdit, QCheckBox, QPushButton, QTableWidget, QTableWidgetItem, QSizePolicy, \
    QTabWidget, QStackedWidget, QTreeView


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

    def checkbox_toggled(self):
        self.data.checked = self.checkbox.isChecked()

    # override default mousePressEvent behaviour
    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)


class SuppFileListItem(UIListItem):
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
        available_width = self.main_window.supp_files_view.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)

    def preview_file(self):
        self.preview_requested.emit(self.file_url)


class View(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        self.multi_page = QStackedWidget(self)
        self.setCentralWidget(self.multi_page)

        self.search_page = QWidget(self)
        self.results_page = QWidget(self)
        
        self.multi_page.addWidget(self.search_page)
        self.multi_page.addWidget(self.results_page) 

        self.init_search_ui_elements()
        self.init_search_layouts()

        self.init_results_ui_elements()
        self.init_results_layouts()

        self.init_load_animation()

    def init_search_ui_elements(self):
        self.search_status = QLabel(self.search_page)
        self.query_field = QLineEdit(self.search_page)
        self.prog_bar = QProgressBar(self.search_page)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()
        self.search_btn = QPushButton("Search", self.search_page)
        self.stop_search_btn = QPushButton("Stop Search", self.search_page)
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", self.search_page)
        self.paper_list = QListWidget(self.search_page)

        self.title_disp = QTextEdit(self.search_page)
        self.title_disp.setPlaceholderText("Title will be shown here")
        self.abstract_disp = QTextEdit(self.search_page)
        self.abstract_disp.setPlaceholderText("Abstract will be shown here")
        self.supp_files_view = QListWidget(self.search_page)
        self.loading_label = QLabel(self.search_page)
        self.loading_label.setAlignment(Qt.AlignCenter)
        self.previews = QWidget(self.search_page)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)

    def init_search_layouts(self):
        pane_0 = QVBoxLayout()
        pane_0.addWidget(QLabel("Enter Query:"))
        pane_0.addWidget(self.query_field)
        pane_0.addWidget(self.prog_bar)
        pane_0.addWidget(self.search_btn)
        pane_0.addWidget(self.stop_search_btn)
        pane_0.addWidget(self.paper_list)
        pane_0.addWidget(self.proceed_btn)

        pane_1 = QVBoxLayout()
        pane_1.addWidget(QLabel("Title:"))
        pane_1.addWidget(self.title_disp)
        pane_1.addWidget(QLabel("Abstract:"))
        pane_1.addWidget(self.abstract_disp)
        pane_1.addWidget(QLabel("Supplementary Files:"))
        pane_1.addWidget(self.supp_files_view)
        pane_1.addWidget(self.loading_label)

        main_pane = QHBoxLayout(self.search_page)
        main_pane.addLayout(pane_0)
        main_pane.addLayout(pane_1)
        main_pane.addWidget(self.previews)
        self.search_page.setLayout(main_pane)

    def init_results_ui_elements(self):
        self.paper_list = QListWidget(self.results_page)
        self.processed_files_view = QListWidget(self.results_page)
        self.previews = QWidget(self.results_page)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)

    def init_results_layouts(self):
        pane_0 = QVBoxLayout()
        pane_0.addWidget(self.paper_list)

        pane_1 = QVBoxLayout()
        pane_1.addWidget(self.processed_files_view)

        main_pane = QHBoxLayout(self.results_page)
        main_pane.addLayout(pane_0)
        main_pane.addLayout(pane_1)
        main_pane.addWidget(self.previews)
        self.results_page.setLayout(main_pane)

    def init_load_animation(self):
        self.load_timer = QTimer(self)
        self.load_dots = 0
        self.load_timer.timeout.connect(self.update_load_text)

    def start_load_animation(self):
        self.load_dots = 0
        self.loading_label.setText("Loading.")
        self.load_timer.start(500)  # ms

    def stop_load_animation(self):
        self.load_timer.stop()
        self.loading_label.clear()

    def update_load_text(self):
        self.load_dots = (self.load_dots + 1) % 4
        self.loading_label.setText("Loading" + "." * self.load_dots)

    def display_paper(self, paper_data, progress):
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

        self.supp_files_view.clear()
        for f in paper.files:
            item_container = QListWidgetItem()
            file_item = SuppFileListItem(f.url, self, f)
            file_item.checkbox.setChecked(f.checked)
            item_container.setSizeHint(file_item.sizeHint())
            self.supp_files_view.addItem(item_container)
            self.supp_files_view.setItemWidget(item_container, file_item)

    def display_multisheet_table(self, df_dict):
        tab_widget = QTabWidget(self.previews)
        for sheet, df in df_dict.items():
            table = self._create_ui_table(df)
            tab_widget.addTab(table, sheet)
        self.previews_layout.addWidget(tab_widget)

    def _create_ui_table(self, data):
        table = QTableWidget()
        table.setColumnCount(len(data.columns))
        table.setHorizontalHeaderLabels(data.columns.astype(str))
        for row, data in data.iterrows():
            table.insertRow(row)
            for col, value in enumerate(data):
                table.setItem(row, col, QTableWidgetItem(str(value)))

        return table

    def switch_to_results_page(self):
        self.multi_page.setCurrentWidget(self.results_page)