from PyQt5.QtCore import pyqtSignal, QRect
from PyQt5.QtWidgets import QStyleOptionButton, QHeaderView, QStyle

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
