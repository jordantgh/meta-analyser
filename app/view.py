from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, pyqtSignal, QTimer, QRect
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QPushButton, QListWidget, QListWidgetItem, QProgressBar, QTextEdit, QCheckBox, QPushButton, QTableWidget, QTableWidgetItem, QSizePolicy,QTabWidget, QStackedWidget, QHeaderView, QStyleOptionButton, QStyle, QSplitter


class CheckableHeaderView(QHeaderView):
    def __init__(self, orientation, parent=None):
        super().__init__(orientation, parent)
        self.setSectionsClickable(True)
        self._checkedSections = set()

    def paintSection(self, painter, rect, logicalIndex):
        painter.save()
        super().paintSection(painter, rect, logicalIndex)
        painter.restore()

        if logicalIndex in self._checkedSections:
            state = QStyle.State_On
        else:
            state = QStyle.State_Off

        option = QStyleOptionButton()
        option.rect = QRect(rect.x() + rect.width() // 2 - 6, rect.y() + rect.height() // 2 - 6, 12, 12)
        option.state = QStyle.State_Enabled | state
        self.style().drawControl(QStyle.CE_CheckBox, option, painter)

    def mousePressEvent(self, event):
        index = self.logicalIndexAt(event.pos())
        if index in self._checkedSections:
            self._checkedSections.remove(index)
        else:
            self._checkedSections.add(index)
        self.viewport().update()
        super().mousePressEvent(event)
        
    def set_all_sections_checked(self):
        self._checkedSections = set(range(self.count()))
        self.viewport().update()


class UIListItem(QWidget):
    def __init__(self, data, title):
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

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        list_widget = self.parent().parent()
        list_item = list_widget.itemAt(self.parent().mapToParent(event.pos()))
        list_widget.setCurrentItem(list_item)


class ArticleListItem(UIListItem):
    def __init__(self, article_data):
        self.article_id = article_data.pmc_id
        super().__init__(article_data, article_data.title)


class DataListItem(UIListItem):
    preview_requested = pyqtSignal(object)

    def __init__(self, main_window, file_data, disp_name_func):
        self.main_window = main_window
        self.file_id = file_data.id
        disp_name = disp_name_func(file_data)
        super().__init__(file_data, disp_name)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        self.preview_requested.emit(self.file_id)

    def get_disp_name(self, text):
        font_metrics = QFontMetrics(self.main_window.font())
        available_width = self.main_window.supp_files_view.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)


class SuppFileListItem(DataListItem):
    def __init__(self, main_window, file_data):
        super().__init__(main_window, file_data, lambda fd: self.get_disp_name(fd.url.split('/')[-1]))
        self.file_url = file_data.url


class ProcessedTableListItem(DataListItem):
    def __init__(self, main_window, file_data):
        super().__init__(main_window, file_data, lambda fd: self.get_disp_name(fd.id))


class View(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        self.multi_page = QStackedWidget(self)
        self.setCentralWidget(self.multi_page)

        self.search_page = QWidget(self)
        self.multi_page.addWidget(self.search_page)

        self.init_search_ui_elements()
        self.init_search_layouts()
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
        self.stop_search_btn.hide()
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", self.search_page)
        self.article_list = QListWidget(self.search_page)

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
 
        self.query_filter_field = QLineEdit(self.search_page)
        self.filter_btn = QPushButton("Filter", self.search_page)

    def init_search_layouts(self):
        pane_0 = QVBoxLayout()
        pane_0.addWidget(QLabel("Enter Query:"))
        pane_0.addWidget(self.query_field)
        pane_0.addWidget(self.prog_bar)
        pane_0.addWidget(self.search_btn)
        pane_0.addWidget(self.stop_search_btn)
        pane_0.addWidget(self.search_status)
        pane_0.addWidget(self.article_list)
        pane_0.addWidget(self.proceed_btn)
        pane_0.addWidget(QLabel("Filter Query:"))
        pane_0.addWidget(self.query_filter_field)
        pane_0.addWidget(self.filter_btn)
        widget_0 = QWidget()
        widget_0.setLayout(pane_0)

        pane_1 = QVBoxLayout()
        pane_1.addWidget(QLabel("Title:"))
        pane_1.addWidget(self.title_disp)
        pane_1.addWidget(QLabel("Abstract:"))
        pane_1.addWidget(self.abstract_disp)
        pane_1.addWidget(QLabel("Supplementary Files:"))
        pane_1.addWidget(self.supp_files_view)
        pane_1.addWidget(self.loading_label)
        pane_1.addWidget(self.prune_btn)
        widget_1 = QWidget()
        widget_1.setLayout(pane_1)

        # Use QSplitter for the main layout
        main_splitter = QSplitter(Qt.Horizontal, self.search_page)
        main_splitter.addWidget(widget_0)
        main_splitter.addWidget(widget_1)
        main_splitter.addWidget(self.previews)

        main_pane = QVBoxLayout(self.search_page)
        main_pane.addWidget(main_splitter)
        self.search_page.setLayout(main_pane)

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

    def suppfilelistitem_factory(self, file_data):
        return SuppFileListItem(self, file_data)
      
    def processedtablelistitem_factory(self, file_data):
        return ProcessedTableListItem(self, file_data)

    def display_article(self, article_data, progress):
        item = QListWidgetItem()
        article_widget = ArticleListItem(article_data)
        item.setSizeHint(article_widget.sizeHint())
        item.setData(Qt.UserRole, article_data.pmc_id)
        self.article_list.addItem(item)
        self.article_list.setItemWidget(item, article_widget)
        self.prog_bar.setValue(progress + 1)

    def clear_article_list_and_files_view(self):
        print("Clearing article list and files view...")
        self.article_list.clear()
        self.supp_files_view.clear()

    def update_article_display(self, article, element_type, list_item_func):
        self.title_disp.setText(article.title)
        self.abstract_disp.setText(article.abstract)

        self.supp_files_view.clear()
        for file_data in getattr(article, element_type):
            item_container = QListWidgetItem()
            file_item = list_item_func(file_data)
            file_item.checkbox.setChecked(file_data.checked)
            item_container.setSizeHint(file_item.sizeHint())
            item_container.setData(Qt.UserRole, file_data.id)
            self.supp_files_view.addItem(item_container)
            self.supp_files_view.setItemWidget(item_container, file_item)

    def populate_filtered_article_list(self, articles, list_item_func):
        print("Populating filtered article list...")
        for article in articles:
            print("Adding article: " + article.pmc_id)
            item = QListWidgetItem()
            article_widget = ArticleListItem(article)
            item.setSizeHint(article_widget.sizeHint())
            item.setData(Qt.UserRole, article.pmc_id)
            self.article_list.addItem(item)
            self.article_list.setItemWidget(item, article_widget)
            self.update_article_display(article, 'processed_tables', list_item_func)

    def display_multisheet_table(self, df_dict, use_checkable_header):
        tab_widget = QTabWidget(self.previews)
        for sheet, df in df_dict.items():
            table = self._create_ui_table(df, use_checkable_header)
            tab_widget.addTab(table, sheet)
        self.previews_layout.addWidget(tab_widget)

    def _create_ui_table(self, data, use_checkable_header):
        table = QTableWidget()
        table.setRowCount(len(data.index))
        table.setColumnCount(len(data.columns))
        table.setHorizontalHeaderLabels(data.columns.astype(str))

        if use_checkable_header:
            header = CheckableHeaderView(Qt.Horizontal, table)
            table.setHorizontalHeader(header)
            header.set_all_sections_checked()
        else:
            header = QHeaderView(Qt.Horizontal, table)
            table.setHorizontalHeader(header)

        for row, row_data in data.iterrows():
            for col, value in enumerate(row_data):
                table.setItem(row, col, QTableWidgetItem(str(value)))

        return table
