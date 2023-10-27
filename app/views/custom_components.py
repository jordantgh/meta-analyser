from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from PyQt5.QtGui import (
        QFocusEvent, QMouseEvent, QStandardItemModel, QPainter
    )

    from PyQt5.QtWidgets import QTabWidget
    from utils.constants import PageIdentity

from PyQt5.QtGui import QKeySequence


class CustomTabBar(QTabBar):
    def focusInEvent(self, event: QFocusEvent):
        if event.reason() == Qt.TabFocusReason:
            event.ignore()
        else:
            super().focusInEvent(event)


class TabPage(QWidget):
    def __init__(self, parent, page_identity):
        super().__init__(parent)
        self.page_identity = page_identity


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
