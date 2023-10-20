from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QLabel, QListWidgetItem, QTableWidget, QTableWidgetItem, QTabWidget, QHeaderView, QSplitter, QAction, QMenu

from views.custom_components import CustomTabBar, TabPage, CheckableHeaderView
from views.list import ArticleListItem, SuppFileListItem, ProcessedTableListItem
from views.page import SearchPageElements, ProcessedPageElements

from utils.constants import PageIdentity


class View(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(1024, 768)
        with open("app/views/styles.qss", "r") as f:
            self.setStyleSheet(f.read())

        self.menu_bar = self.menuBar()
        self.file_menu = QMenu("File", self)
        self.save_action = QAction("Save", self)
        self.load_action = QAction("Load", self)
        self.file_menu.addAction(self.save_action)
        self.file_menu.addAction(self.load_action)
        self.menu_bar.addMenu(self.file_menu)

        self.tab_widget = QTabWidget(self)
        self.tab_widget.setTabBar(CustomTabBar())
        self.setCentralWidget(self.tab_widget)

        self.search_tab = TabPage(self, PageIdentity.SEARCH)
        self.parsed_tab = TabPage(self, PageIdentity.PARSED)
        self.pruned_tab = TabPage(self, PageIdentity.PRUNED)

        self.tab_widget.addTab(self.search_tab, "Search")
        self.tab_widget.addTab(self.parsed_tab, "Parsing Results")
        self.tab_widget.addTab(self.pruned_tab, "Pruned Results")

        self.search_elems = SearchPageElements(self.search_tab)
        self.parsed_elems = ProcessedPageElements(self.parsed_tab)
        self.pruned_elems = ProcessedPageElements(self.pruned_tab)

        self.init_search_layouts(self.search_elements)
        self.init_processed_page_layouts(self.parsed_tab, self.parsed_elements)
        self.init_processed_page_layouts(self.pruned_tab, self.pruned_elements)

        self.search_elements.query_field.setFocus()
        self.init_load_animation()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Escape:
            self.focusWidget().clearFocus()

    # getters for active page
    @property
    def active_tab(self):
        return self.tab_widget.currentWidget()

    @property
    def active_elements(self):
        return {
            self.search_tab: self.search_elements,
            self.parsed_tab: self.parsed_elements,
            self.pruned_tab: self.pruned_elements
        }.get(self.active_tab)

    # set active page
    def set_active_tab(self, page_identity):
        if page_identity == PageIdentity.SEARCH:
            self.tab_widget.setCurrentWidget(self.search_tab)
        elif page_identity == PageIdentity.PARSED:
            self.tab_widget.setCurrentWidget(self.parsed_tab)
        elif page_identity == PageIdentity.PRUNED:
            self.tab_widget.setCurrentWidget(self.pruned_tab)

    def init_search_layouts(self, elements):
        left_pane = QVBoxLayout()
        left_pane.addWidget(QLabel("Enter Query:"))
        left_pane.addWidget(elements.query_field)
        left_pane.addWidget(elements.prog_bar)
        left_pane.addWidget(elements.search_btn)
        left_pane.addWidget(elements.stop_search_btn)
        left_pane.addWidget(elements.search_status)
        left_pane.addWidget(elements.article_list_view)
        left_pane.addWidget(QLabel("Associated Data:"))
        left_pane.addWidget(elements.data_list_view)
        left_pane.addWidget(elements.proceed_btn)

        left_pane.setStretchFactor(elements.article_list_view, 3)
        left_pane.setStretchFactor(elements.data_list_view, 1)

        self.init_core_layouts(self.search_tab, elements, left_pane)

    def init_processed_page_layouts(self, page, elements):
        left_pane = QVBoxLayout()
        left_pane.addWidget(elements.prog_bar)
        left_pane.addWidget(elements.article_list_view)
        left_pane.addWidget(QLabel("Associated Data:"))
        left_pane.addWidget(elements.data_list_view)
        left_pane.addWidget(QLabel("Filter Query:"))
        left_pane.addWidget(elements.query_filter_field)
        left_pane.addWidget(elements.filter_btn)
        left_pane.addWidget(elements.prune_btn)

        left_pane.setStretchFactor(elements.article_list_view, 3)
        left_pane.setStretchFactor(elements.data_list_view, 1)

        self.init_core_layouts(page, elements, left_pane)

    def init_core_layouts(self, page, elements, left_pane):
        widget_0 = QWidget()
        widget_0.setLayout(left_pane)

        mid_pane = QVBoxLayout()
        textbox_label = QLabel("Title/Abstract:")
        mid_pane.addWidget(textbox_label)
        mid_pane.addWidget(elements.title_abstract_disp)
        mid_pane.setStretchFactor(textbox_label, 0)
        mid_pane.setStretchFactor(elements.title_abstract_disp, 1)
        mid_widget = QWidget()
        mid_widget.setLayout(mid_pane)

        # Create a container for previews
        preview_pane = QVBoxLayout()
        preview_label = QLabel("Data Preview:")
        preview_pane.addWidget(preview_label)
        preview_pane.addWidget(elements.previews)
        preview_pane.addWidget(elements.loading_label)
        preview_pane.setStretchFactor(preview_label, 0)
        preview_pane.setStretchFactor(elements.previews, 1)
        preview_pane.setStretchFactor(elements.loading_label, 0)
        preview_widget = QWidget()
        preview_widget.setLayout(preview_pane)

        # Create a QSplitter for the title/abstract and previews
        mid_splitter = QSplitter(Qt.Vertical)
        mid_splitter.addWidget(mid_widget)
        mid_splitter.addWidget(preview_widget)
        mid_splitter.setStretchFactor(1,10)


        main_splitter = QSplitter(Qt.Horizontal, page)
        main_splitter.addWidget(widget_0)
        main_splitter.addWidget(mid_splitter)
        main_splitter.setSizes([400, 600])

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
            "LOADING" + "." * self.load_dots
        )

    def show_searching_view(self):
        self.search_elements.prog_bar.setValue(0)
        self.search_elements.prog_bar.show()
        self.search_elements.search_status.setText("Searching...")
        self.search_elements.search_status.show()
        self.search_elements.stop_search_btn.show()
        self.search_elements.stop_search_btn.setEnabled(True)
        
    def hide_searching_view(self):
        self.search_elements.search_status.setText("Stopping search...")
        self.search_elements.prog_bar.hide()
        self.search_elements.search_status.setText("Search stopped.")
        self.search_elements.stop_search_btn.hide()
        self.search_elements.stop_search_btn.setEnabled(False)

    def to_list(self, list_widget, item_widget, id):
        item = QListWidgetItem()
        item.setSizeHint(item_widget.sizeHint())
        item.setData(Qt.UserRole, id)
        list_widget.addItem(item)
        list_widget.setItemWidget(item, item_widget)

    def display_article(self, elements, article, progress):
        article_item = ArticleListItem(article, elements.page_identity)
        self.to_list(elements.article_list_view, article_item, article.pmc_id)
        elements.prog_bar.setValue(progress + 1)

    def update_article_display(self, article, elements, data_set):
        self.clear_list_and_observers(elements.data_list_view)
        elements.title_abstract_disp.setHtml(
            f"<a href='{article.url}'>"
            f"<b>{article.title}</b></a>"
            f"<br><br>{article.abstract}"
        )

        for data in data_set:
            data_item = self.list_item_factory(data, elements.page_identity)
            data_item.checkbox.setChecked(data.checked)
            self.to_list(elements.data_list_view, data_item, data.id)

    def clear_page_lists(self, elements):
        self.clear_list_and_observers(elements.article_list_view)
        self.clear_list_and_observers(elements.data_list_view)

    def clear_list_and_observers(self, list_widget):
        for index in range(list_widget.count()):
            item = list_widget.item(index)
            if item:
                widget = list_widget.itemWidget(item)
                if widget and widget.data.alert_observers():
                    widget.remove()

        list_widget.clear()

    # fully reset the view
    def reset(self):
        for elems in [self.search_tab, self.parsed_tab, self.pruned_tab]:
            self.clear_page_lists(elems)
            elems.title_abstract_disp.clear()
            elems.previews.clear()

    def list_item_factory(self, file_data, context):
        if context == PageIdentity.SEARCH:
            return SuppFileListItem(self, file_data, context)
        else:
            return ProcessedTableListItem(self, file_data, context)

    def display_multisheet_table(
        self,
        df_dict,
        use_checkable_header,
        table_id=None,
        callback=None,
        checked_columns=None
    ):
        tab_widget = self.active_elements.previews
        tab_widget.clear()

        for sheet, df in df_dict.items():
            table = self._create_ui_table(
                df, use_checkable_header, table_id, callback, checked_columns
            )
            tab_widget.addTab(table, sheet)

    def _create_ui_table(
        self,
        data,
        use_checkable_header,
        table_id=None,
        callback=None,
        checked_columns=None
    ):
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
