from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, pyqtSignal, QTimer, QRect
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLineEdit, QLabel, QPushButton, QListWidget, QListWidgetItem, QProgressBar, QTextEdit, QCheckBox, QPushButton, QTableWidget, QTableWidgetItem, QSizePolicy,QTabWidget, QStackedWidget, QHeaderView, QStyleOptionButton, QStyle, QSplitter


class CheckableHeaderView(QHeaderView):
    columns_checked = pyqtSignal(object, list)

    def __init__(self, orientation, parent=None):
        super().__init__(orientation, parent)
        self.setSectionsClickable(True)
        self._checkedSections = set()
        self.table_id = None

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
        self.columns_checked.emit(self.table_id, list(self._checkedSections))
        super().mousePressEvent(event)
        
    def set_checked_sections(self, checked_sections):
        self._checkedSections = set(checked_sections)
        self.viewport().update()
    
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
        self.page = main_window.active_elements
        self.file_id = file_data.id
        disp_name = disp_name_func(file_data)
        super().__init__(file_data, disp_name)

    def mousePressEvent(self, event):
        super().mousePressEvent(event)
        self.preview_requested.emit(self.file_id)

    def get_disp_name(self, text):
        font_metrics = QFontMetrics(self.main_window.font())
        available_width = self.page.supp_files_view.width() - 150
        return font_metrics.elidedText(text, Qt.ElideMiddle, available_width)


class SuppFileListItem(DataListItem):
    def __init__(self, main_window, file_data):
        super().__init__(main_window, file_data, lambda fd: self.get_disp_name(fd.url.split('/')[-1]))
        self.file_url = file_data.url


class ProcessedTableListItem(DataListItem):
    def __init__(self, main_window, file_data):
        super().__init__(main_window, file_data, lambda fd: self.get_disp_name(fd.id))


class CommonPageElements:
    def __init__(self, page):
        self.prog_bar = QProgressBar(page)
        self.prog_bar.setRange(0, 100)
        self.prog_bar.setValue(0)
        self.prog_bar.hide()
        self.article_list = QListWidget(page)
        self.title_disp = QTextEdit(page)
        self.title_disp.setPlaceholderText("Title will be shown here")
        self.abstract_disp = QTextEdit(page)
        self.abstract_disp.setPlaceholderText("Abstract will be shown here")
        self.supp_files_view = QListWidget(page)
        self.previews = QWidget(page)
        self.previews_layout = QVBoxLayout(self.previews)
        self.previews.setLayout(self.previews_layout)
        self.loading_label = QLabel(page)
        self.loading_label.setAlignment(Qt.AlignCenter)

class SearchPageElements(CommonPageElements):
    def __init__(self, page):
        super().__init__(page)
        self.search_status = QLabel(page)
        self.query_field = QLineEdit(page)
        self.search_btn = QPushButton("Search", page)
        self.stop_search_btn = QPushButton("Stop Search", page)
        self.stop_search_btn.hide()
        self.stop_search_btn.setEnabled(False)
        self.proceed_btn = QPushButton("Proceed", page)
        

class ProcessedPageElements(CommonPageElements):
    def __init__(self, page):
        super().__init__(page)
        self.query_filter_field = QLineEdit(page)
        self.filter_btn = QPushButton("Filter", page)
        self.prune_btn = QPushButton("Prune Tables and Columns", page)        

class View(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        
        self.tab_widget = QTabWidget(self)
        self.setCentralWidget(self.tab_widget)

        self.search_page = QWidget(self)
        self.parsed_page = QWidget(self)
        self.pruned_page = QWidget(self)
        
        self.tab_widget.addTab(self.search_page, "Search")
        self.tab_widget.addTab(self.parsed_page, "Parsing Results")
        self.tab_widget.addTab(self.pruned_page, "Pruned Results")
        
        self.search_components = SearchPageElements(self.search_page)
        self.parsed_components = ProcessedPageElements(self.parsed_page)
        self.pruned_components = ProcessedPageElements(self.pruned_page)

        self.init_search_layouts(self.search_components)
        self.init_processed_page_layouts(self.parsed_page, self.parsed_components)
        self.init_processed_page_layouts(self.pruned_page, self.pruned_components)
        self.init_load_animation()

    @property
    def active_page(self):
        return self.tab_widget.currentWidget()

    @property
    def active_elements(self):
        return {
            self.search_page: self.search_components,
            self.parsed_page: self.parsed_components,
            self.pruned_page: self.pruned_components
        }.get(self.active_page, None)

    def init_search_layouts(self, components):
        left_pane = QVBoxLayout()
        left_pane.addWidget(QLabel("Enter Query:"))
        left_pane.addWidget(components.query_field)
        left_pane.addWidget(components.prog_bar)
        left_pane.addWidget(components.search_btn)
        left_pane.addWidget(components.stop_search_btn)
        left_pane.addWidget(components.search_status)
        left_pane.addWidget(components.article_list)
        left_pane.addWidget(components.proceed_btn)
        
        self.init_core_layouts(self.search_page, components, left_pane)

    def init_processed_page_layouts(self, page, components):
        left_pane = QVBoxLayout()
        left_pane.addWidget(components.prog_bar)
        left_pane.addWidget(components.article_list)
        left_pane.addWidget(QLabel("Filter Query:"))
        left_pane.addWidget(components.query_filter_field)
        left_pane.addWidget(components.filter_btn)       
        left_pane.addWidget(components.prune_btn)

        self.init_core_layouts(page, components, left_pane)
   
    def init_core_layouts(self, page, components, left_pane):        
        widget_0 = QWidget()
        widget_0.setLayout(left_pane)

        mid_pane = QVBoxLayout()
        mid_pane.addWidget(QLabel("Title:"))
        mid_pane.addWidget(components.title_disp)
        mid_pane.addWidget(QLabel("Abstract:"))
        mid_pane.addWidget(components.abstract_disp)
        mid_pane.addWidget(QLabel("Supplementary Files:"))
        mid_pane.addWidget(components.supp_files_view)
        
        widget_1 = QWidget()
        widget_1.setLayout(mid_pane)

        main_splitter = QSplitter(Qt.Horizontal, page)
        main_splitter.addWidget(widget_0)
        main_splitter.addWidget(widget_1)
        main_splitter.addWidget(components.previews)

        main_pane = QVBoxLayout(page)
        main_pane.addWidget(main_splitter)
        page.setLayout(main_pane)

    def init_load_animation(self):
        self.load_timer = QTimer(self)
        self.load_dots = 0
        self.load_timer.timeout.connect(self.update_load_text)

    def start_load_animation(self):
        self.load_dots = 0
        self.active_elements.loading_label.setText("Loading.")
        self.load_timer.start(500)  # ms

    def stop_load_animation(self):
        self.load_timer.stop()
        self.active_elements.loading_label.clear()

    def update_load_text(self):
        self.load_dots = (self.load_dots + 1) % 4
        self.active_elements.loading_label.setText("Loading" + "." * self.load_dots)

    def suppfilelistitem_factory(self, file_data):
        return SuppFileListItem(self, file_data)
      
    def processedtablelistitem_factory(self, file_data):
        return ProcessedTableListItem(self, file_data)

    def display_article(self, article_data, progress):
        item = QListWidgetItem()
        article_widget = ArticleListItem(article_data)
        item.setSizeHint(article_widget.sizeHint())
        item.setData(Qt.UserRole, article_data.pmc_id)
        self.active_elements.article_list.addItem(item)
        self.active_elements.article_list.setItemWidget(item, article_widget)
        self.active_elements.prog_bar.setValue(progress + 1)

    def clear_article_list_and_files_view(self):
        self.active_elements.article_list.clear()
        self.active_elements.supp_files_view.clear()

    def update_article_display(self, article, element_type, list_item_func):
        self.active_elements.supp_files_view.clear()
        
        # This monster must be slain
        if element_type == 'to_prune' and article.to_prune:
            tables_to_process = [table for table in article.processed_tables if table.id in article.to_prune and table.checked]
        elif element_type == 'supp_files':
            tables_to_process = getattr(article, element_type)
        else:
            tables_to_process = [table for table in article.processed_tables if table.checked]

        for data in tables_to_process:
            item_container = QListWidgetItem()
            file_item = list_item_func(data)
            file_item.checkbox.setChecked(data.checked)
            item_container.setSizeHint(file_item.sizeHint())
            item_container.setData(Qt.UserRole, data.id)
            self.active_elements.supp_files_view.addItem(item_container)
            self.active_elements.supp_files_view.setItemWidget(item_container, file_item)

    def display_multisheet_table(self, df_dict, use_checkable_header, table_id=None, callback=None, checked_columns=None):
        tab_widget = QTabWidget(self.active_elements.previews)
        for sheet, df in df_dict.items():
            table = self._create_ui_table(df, use_checkable_header, table_id, callback, checked_columns)
            tab_widget.addTab(table, sheet)
        self.active_elements.previews_layout.addWidget(tab_widget)

    def _create_ui_table(self, data, use_checkable_header, table_id=None, callback=None, checked_columns=None):
        ui_table = QTableWidget()
        ui_table.setRowCount(len(data.index))
        ui_table.setColumnCount(len(data.columns))
        ui_table.setHorizontalHeaderLabels(data.columns.astype(str))

        if use_checkable_header:
            header = CheckableHeaderView(Qt.Horizontal, ui_table)
            header.table_id = table_id
            header.columns_checked.connect(callback)
            ui_table.setHorizontalHeader(header)
            if checked_columns is not None:
                header.set_checked_sections(checked_columns)
            else:
                header.set_all_sections_checked()
        else:
            header = QHeaderView(Qt.Horizontal, ui_table)
            ui_table.setHorizontalHeader(header)

        for row, row_data in data.iterrows():
            for col, value in enumerate(row_data):
                ui_table.setItem(row, col, QTableWidgetItem(str(value)))

        return ui_table
