from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyQt5.QtGui import (
        QFocusEvent, QMouseEvent, QStandardItemModel, QPainter
    )

    from PyQt5.QtWidgets import QTabWidget
    from utils.constants import PageIdentity

from PyQt5.QtGui import QKeySequence
from PyQt5.QtCore import Qt, QEvent, pyqtSignal, QRect
from PyQt5.QtWidgets import (
    QStyleOptionButton, QHeaderView, QStyle, QTabBar, QWidget, QVBoxLayout,
    QTableView, QShortcut, QLineEdit, QDialog, QPushButton, QMessageBox
)


class CustomTabBar(QTabBar):
    def focusInEvent(self, event: 'QFocusEvent'):
        if event.reason() == Qt.TabFocusReason:
            event.ignore()
        else:
            super().focusInEvent(event)


class TabPage(QWidget):
    def __init__(self, parent: 'QTabWidget', page_identity: 'PageIdentity'):
        super().__init__(parent)
        self.page_identity = page_identity


class CustomTable(QTableView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.shortcut = QShortcut(QKeySequence("Ctrl+F"), self)
        self.shortcut.activated.connect(self.showFindDialog)

    def showFindDialog(self):
        findDialog = QDialog(self)
        findDialog.setWindowTitle("Find")
        global_position = self.mapToGlobal(self.pos())
        findDialog.setGeometry(global_position.x() + 50,
                               global_position.y() + 50, 300, 100)

        layout = QVBoxLayout()
        findText = QLineEdit()
        findButton = QPushButton("Find")

        layout.addWidget(findText)
        layout.addWidget(findButton)

        findDialog.setLayout(layout)

        findButton.clicked.connect(lambda: self.findText(findText.text()))

        # Close on 'Esc' key
        findDialog.installEventFilter(self)
        findDialog.exec_()

    def eventFilter(self, obj: 'QDialog', event: 'QEvent') -> 'bool':
        if event.type() == QEvent.KeyRelease and event.key() == Qt.Key_Escape:
            obj.close()
            return True
        return super().eventFilter(obj, event)

    def findText(self, text: 'str'):
        text = text.lower()
        model: 'QStandardItemModel' = self.model()
        for row in range(model.rowCount()):
            for col in range(model.columnCount()):
                item = model.item(row, col)
                if item and item.text().lower() == text:
                    self.selectRow(row)
                    return
        QMessageBox.information(self, "No Matches", "No matches found.")


class CheckableHeaderView(QHeaderView):
    columns_checked = pyqtSignal(object, list)

    def __init__(
        self,
        orientation: 'Qt.Orientation',
        parent: 'QTableView' = None,
        table_id: 'str' = None
    ):
        super().__init__(orientation, parent)
        self.setSectionsClickable(True)
        self._checkedSections = set()
        self.table_id = table_id

    def paintSection(
        self,
        painter: 'QPainter',
        rect: 'QRect',
        logicalIndex: 'int'
    ):
        painter.save()
        super().paintSection(painter, rect, logicalIndex)
        painter.restore()

        if logicalIndex in self._checkedSections:
            state = QStyle.State_On
        else:
            state = QStyle.State_Off

        option = QStyleOptionButton()
        option.rect = QRect(
            rect.x() + rect.width() // 2 - 6,
            rect.y() + rect.height() // 2 - 6, 12, 12
        )

        option.state = QStyle.State_Enabled | state
        self.style().drawControl(QStyle.CE_CheckBox, option, painter)

    def mousePressEvent(self, event: 'QMouseEvent'):
        index = self.logicalIndexAt(event.pos())

        clickable_area = 50

        # Calculate the offsets to center the clickable area
        x_offset = (self.sectionSize(index) - clickable_area) // 2
        y_offset = (self.height() - clickable_area) // 2
        checkbox_rect = QRect(
            self.sectionViewportPosition(index) + x_offset,
            y_offset,
            clickable_area,
            clickable_area
        )

        if checkbox_rect.contains(event.pos()):
            if index in self._checkedSections:
                self._checkedSections.remove(index)
            else:
                self._checkedSections.add(index)
            self.viewport().update()
            self.columns_checked.emit(
                self.table_id, list(self._checkedSections))
        else:
            super().mousePressEvent(event)

    def set_checked_sections(self, checked_sections: 'list[int]'):
        self._checkedSections = set(checked_sections)
        self.viewport().update()

    def set_all_sections_checked(self):
        self._checkedSections = set(range(self.count()))
        self.viewport().update()
