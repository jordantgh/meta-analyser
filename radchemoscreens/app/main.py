from PyQt5.QtWidgets import QApplication
from model import CRISPRModel
from view import CRISPRView
from controller import CRISPRController

def main():
    app = QApplication([])

    model = CRISPRModel()
    view = CRISPRView()
    controller = CRISPRController(model, view)

    view.show()
    app.exec_()

if __name__ == "__main__":
    main()
