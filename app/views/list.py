from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QWidget, QHBoxLayout, QLabel, QCheckBox, QSizePolicy
from PyQt5.QtGui import QFontMetrics

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
        if self.data.alert_observers():
            self.data.checkbox_toggled()

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

        file_data.register_observer(self)

    def remove(self):
        self.data.remove_observer(self)

    def update(self, processed_table):
        # Disconnect the signal before updating the checkbox to avoid triggering the signal again, then reconnect it
        self.checkbox.toggled.disconnect(self.checkbox_toggled)
        self.checkbox.setChecked(processed_table.checked)
        self.checkbox.toggled.connect(self.checkbox_toggled)
