from typing import TYPE_CHECKING, Callable, cast
if TYPE_CHECKING:
    from views.page import PageElements
    from PyQt5.QtGui import QKeyEvent
    from PyQt5.QtWidgets import QListWidget
    from model.article_managers import Article, BaseData
    from views.list import ListItem, DataListItem
    from pandas import DataFrame

from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtCore import Qt, QTimer, QCoreApplication
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QLabel, QListWidgetItem, QTabWidget,
    QHeaderView, QSplitter, QAction, QMenu, QAbstractItemView
)

from views.custom_components import (
    CustomTabBar, TabPage, CustomTable, CheckableHeaderView
)
from views.list import (
    ListItem, ArticleListItem, SuppFileListItem, ProcessedTableListItem
)
from views.page import SearchPageElements, ProcessedPageElements

from utils.constants import PageIdentity

import pandas as pd


class View(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(1024, 768)
        with open("app/views/styles.qss", "r") as f:
            self.setStyleSheet(f.read())

        self.menu_bar = self.menuBar()
        self.file_menu = QMenu("File", self)
        self.save_action = QAction("Save", self)
        self.save_as_action = QAction("Save as...", self)
        self.load_action = QAction("Load", self)
        self.file_menu.addAction(self.save_action)
        self.file_menu.addAction(self.save_as_action)
        self.file_menu.addAction(self.load_action)
        self.save_action.setShortcut("Ctrl+S")
        self.save_as_action.setShortcut("Ctrl+Shift+S")
        self.load_action.setShortcut("Ctrl+O")
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

        self._init_search_layouts(self.search_elems)
        self._init_processed_page_layouts(self.parsed_tab, self.parsed_elems)
        self._init_processed_page_layouts(self.pruned_tab, self.pruned_elems)

        self.search_elems.query_field.setFocus()
        self._init_load_animation()

    def keyPressEvent(self, event: 'QKeyEvent'):
        if event.key() == Qt.Key_Escape:
            self.focusWidget().clearFocus()

    # Getters for active page
    @property
    def active_tab(self) -> 'TabPage':
        return self.tab_widget.currentWidget()

    @property
    def active_elements(self) -> 'PageElements':
        return {
            self.search_tab: self.search_elems,
            self.parsed_tab: self.parsed_elems,
            self.pruned_tab: self.pruned_elems
        }.get(self.active_tab)

    # Set active page
    def set_active_tab(self, page_identity: 'PageIdentity'):
        if page_identity == PageIdentity.SEARCH:
            self.tab_widget.setCurrentWidget(self.search_tab)
        elif page_identity == PageIdentity.PARSED:
            self.tab_widget.setCurrentWidget(self.parsed_tab)
        elif page_identity == PageIdentity.PRUNED:
            self.tab_widget.setCurrentWidget(self.pruned_tab)

    def _init_search_layouts(self, elements: 'SearchPageElements'):
        left_pane = QVBoxLayout()
        left_pane.addWidget(QLabel("Enter Query:"))
        left_pane.addWidget(elements.query_field)
        left_pane.addWidget(elements.prog_bar)
        left_pane.addWidget(elements.search_btn)
        left_pane.addWidget(elements.stop_search_btn)
        left_pane.addWidget(elements.search_status)
        left_pane.addWidget(elements.article_ui_list)
        left_pane.addWidget(QLabel("Associated Data:"))
        left_pane.addWidget(elements.data_ui_list)
        left_pane.addWidget(elements.proceed_btn)

        left_pane.setStretchFactor(elements.article_ui_list, 3)
        left_pane.setStretchFactor(elements.data_ui_list, 1)

        self._init_core_layouts(self.search_tab, elements, left_pane)

    def _init_processed_page_layouts(
            self,
            page: 'TabPage',
            elements: 'ProcessedPageElements'
    ):
        left_pane = QVBoxLayout()
        left_pane.addWidget(elements.prog_bar)
        left_pane.addWidget(elements.article_ui_list)
        left_pane.addWidget(QLabel("Associated Data:"))
        left_pane.addWidget(elements.data_ui_list)
        left_pane.addWidget(QLabel("Filter Query:"))
        left_pane.addWidget(elements.query_filter_field)
        left_pane.addWidget(elements.filter_btn)
        left_pane.addWidget(elements.prune_btn)

        left_pane.setStretchFactor(elements.article_ui_list, 3)
        left_pane.setStretchFactor(elements.data_ui_list, 1)

        self._init_core_layouts(page, elements, left_pane)

    def _init_core_layouts(
        self,
        page: 'TabPage',
        elements: 'PageElements',
        left_pane: 'QVBoxLayout'
    ):
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

        # Container for previews
        preview_pane = QVBoxLayout()
        preview_pane.addWidget(elements.outer_tab_widget)
        preview_pane.addWidget(elements.loading_label)
        preview_pane.setStretchFactor(elements.outer_tab_widget, 1)
        preview_pane.setStretchFactor(elements.loading_label, 0)
        preview_widget = QWidget()
        preview_widget.setLayout(preview_pane)

        # QSplitter for the title/abstract and previews
        mid_splitter = QSplitter(Qt.Vertical)
        mid_splitter.addWidget(mid_widget)
        mid_splitter.addWidget(preview_widget)
        mid_splitter.setStretchFactor(1, 10)

        main_splitter = QSplitter(Qt.Horizontal, page)
        main_splitter.addWidget(widget_0)
        main_splitter.addWidget(mid_splitter)
        main_splitter.setSizes([400, 600])

        main_pane = QVBoxLayout(page)
        main_pane.addWidget(main_splitter)
        page.setLayout(main_pane)

    def _init_load_animation(self):
        self.load_timer = QTimer(self)
        self.load_dots = 0
        self.load_timer.timeout.connect(self._update_load_text)

    def start_load_animation(self):
        self.load_dots = 0
        self.active_elements.loading_label.setText("Loading.")
        self.load_timer.start(500)  # ms

    def stop_load_animation(self):
        self.load_timer.stop()
        self.active_elements.loading_label.clear()

    def _update_load_text(self):
        self.load_dots = (self.load_dots + 1) % 4
        self.active_elements.loading_label.setText(
            "LOADING" + "." * self.load_dots
        )

    def show_searching_view(self):
        self.search_elems.prog_bar.setValue(0)
        self.search_elems.prog_bar.show()
        self.search_elems.search_status.setText("Searching...")
        self.search_elems.search_status.show()
        self.search_elems.stop_search_btn.show()
        self.search_elems.stop_search_btn.setEnabled(True)

    def hide_searching_view(self):
        self.search_elems.search_status.setText("Stopping search...")
        self.search_elems.prog_bar.hide()
        self.search_elems.search_status.setText("Search stopped.")
        self.search_elems.stop_search_btn.hide()
        self.search_elems.stop_search_btn.setEnabled(False)

    def _to_list(
        self,
        list: 'QListWidget',
        item_widget: 'QWidget',
        id: 'str'
    ):
        item = QListWidgetItem()
        item.setSizeHint(item_widget.sizeHint())
        item.setData(Qt.UserRole, id)
        list.addItem(item)
        list.setItemWidget(item, item_widget)

    def display_article(
        self,
        elements: 'PageElements',
        article: 'Article',
        progress: 'int'
    ):
        article_item = ArticleListItem(article, elements.page_identity)
        self._to_list(elements.article_ui_list, article_item, article.pmc_id)
        elements.prog_bar.setValue(progress + 1)

    def update_article_display(
        self,
        article: 'Article',
        elements: 'PageElements',
        data_set: 'list[BaseData]'
    ):
        self.clear_list_and_observers(elements.data_ui_list)
        elements.title_abstract_disp.setHtml(
            f"<a href='{article.url}'>"
            f"<b>{article.title}</b></a>"
            f"<br><br>{article.abstract}"
        )

        for data in data_set:
            data_item: 'DataListItem' = self._list_item_factory(
                data, elements.page_identity
            )
            data_item.checkbox.setChecked(data.checked)
            self._to_list(elements.data_ui_list, data_item, data.id)

    def clear_page_lists(self, elements: 'PageElements'):
        self.clear_list_and_observers(elements.article_ui_list)
        self.clear_list_and_observers(elements.data_ui_list)

    def clear_list_and_observers(self, list_widget: 'QListWidget'):
        for index in range(list_widget.count()):
            item = list_widget.item(index)
            if item:
                widget = cast(ListItem, list_widget.itemWidget(item))
                if widget and widget.data.alert_observers():
                    widget.remove()

        list_widget.clear()

    def reset(self):
        elems: 'PageElements'
        for elems in [self.search_elems, self.parsed_elems, self.pruned_elems]:
            self.clear_page_lists(elems)
            elems.title_abstract_disp.clear()
            elems.previews.clear()

    def _list_item_factory(
        self, file_data: 'BaseData', context: 'PageIdentity'
    ) -> 'DataListItem':
        if context == PageIdentity.SEARCH:
            return SuppFileListItem(self, file_data, context)
        else:
            return ProcessedTableListItem(self, file_data, context)

    def display_multisheet_table(
        self,
        df_dict: 'dict[str, DataFrame]',
        use_checkable_header: 'bool',
        table_id: 'str' = None,
        callback: 'Callable' = None,
        checked_columns: 'list[int]' = None
    ):
        tab_widget = self.active_elements.previews
        tab_widget.clear()

        for sheet, df in df_dict.items():
            table = self._create_ui_table(
                df,
                use_checkable_header,
                table_id,
                callback,
                checked_columns
            )
            tab_widget.addTab(table, sheet)

    def _create_ui_table(
        self,
        data: 'DataFrame',
        use_checkable_header: 'bool',
        table_id: 'str' = None,
        callback: 'Callable' = None,
        checked_columns: 'list[int]' = None
    ) -> 'CustomTable':

        model = QStandardItemModel()
        model.setHorizontalHeaderLabels(data.columns.astype(str))

        # Disable GUI updates
        ui_table = CustomTable()
        ui_table.setUpdatesEnabled(False)

        # Populate the table in chunks
        chunk_size = 5000
        for i in range(0, len(data), chunk_size):
            chunk = data.iloc[i:i + chunk_size]
            for row in chunk.itertuples(index=False, name=None):
                items = [
                    QStandardItem(str(value)) if not pd.isna(value) else
                    QStandardItem("") for value in row
                ]
                model.appendRow(items)
            QCoreApplication.processEvents()

        ui_table.setModel(model)

        # Re-enable GUI updates
        ui_table.setUpdatesEnabled(True)
        ui_table.update()

        # Set Read-Only
        ui_table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        # Disable Column Resizing
        ui_table.horizontalHeader().setSectionResizeMode(QHeaderView.Fixed)

        # Checkable Header
        if use_checkable_header:
            header = CheckableHeaderView(Qt.Horizontal, ui_table, table_id)
            header.columns_checked.connect(callback)
            ui_table.setHorizontalHeader(header)
            if checked_columns is not None:
                header.set_checked_sections(checked_columns)
            else:
                header.set_all_sections_checked()

        return ui_table
