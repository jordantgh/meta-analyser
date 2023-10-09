from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QLabel, QListWidgetItem, QTableWidget, QTableWidgetItem, QTabWidget, QHeaderView, QSplitter

from views.custom_components import CheckableHeaderView
from views.list import ArticleListItem, SuppFileListItem, ProcessedTableListItem
from views.page import SearchPageElements, ProcessedPageElements


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
        self.init_processed_page_layouts(
            self.parsed_page, self.parsed_components)
        self.init_processed_page_layouts(
            self.pruned_page, self.pruned_components)
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
        self.active_elements.loading_label.setText(
            "Loading" + "." * self.load_dots)

    # the context argument is only used as a dummy in this specific function
    # because the it's used in the update_article_display function later
    # TODO is to harmonise this so we arent passing around dummy arguments
    def suppfilelistitem_factory(self, file_data, context):
        return SuppFileListItem(self, file_data, context)

    def processedtablelistitem_factory(self, file_data, context):
        return ProcessedTableListItem(self, file_data, context)

    def display_article(self, components, context, article_data, progress):
        item = QListWidgetItem()
        article_widget = ArticleListItem(article_data, context)
        item.setSizeHint(article_widget.sizeHint())
        item.setData(Qt.UserRole, article_data.pmc_id)
        components.article_list.addItem(item)
        components.article_list.setItemWidget(item, article_widget)
        components.prog_bar.setValue(progress + 1)

    def clear_supp_files_view(self):
        # here we get all the list items and deregister their observers
        # removing them
        for i in range(self.active_elements.supp_files_view.count()):
            item = self.active_elements.supp_files_view.item(i)

            if item:
                list_item = self.active_elements.supp_files_view.itemWidget(
                    item)

                if list_item and list_item.data.alert_observers():
                    list_item.remove()

        self.active_elements.supp_files_view.clear()

    def clear_article_list_and_files_view(self):
        self.active_elements.article_list.clear()
        self.clear_supp_files_view()

    def update_article_display(self, article, element_type, list_item_func, context):
        self.clear_supp_files_view()
        self.active_elements.title_disp.setText(article.title)
        self.active_elements.abstract_disp.setText(article.abstract)

        # TODO this monster must be slain
        tables_to_display = []
        if element_type == 'pruned_article_tables':
            tables_to_display = article.pruned_tables
        elif element_type == 'supp_files':
            tables_to_display = article.supp_files
        else:
            tables_to_display = article.processed_tables

        for data in tables_to_display:
            item_container = QListWidgetItem()
            file_item = list_item_func(data, context)
            file_item.checkbox.setChecked(data.checked)
            item_container.setSizeHint(file_item.sizeHint())
            item_container.setData(Qt.UserRole, data.id)
            self.active_elements.supp_files_view.addItem(item_container)
            self.active_elements.supp_files_view.setItemWidget(
                item_container, file_item)

    def display_multisheet_table(self, df_dict, use_checkable_header, table_id=None, callback=None, checked_columns=None):
        tab_widget = QTabWidget(self.active_elements.previews)
        for sheet, df in df_dict.items():
            table = self._create_ui_table(
                df, use_checkable_header, table_id, callback, checked_columns)
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
