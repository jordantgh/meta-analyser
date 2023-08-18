from PyQt5.QtWidgets import QApplication
from model import Model
from view import View
from controller import Controller

def main():
    app = QApplication([])

    model = Model()
    view = View()
    _ = Controller(model, view)

    view.show()
    app.exec_()

if __name__ == "__main__":
    main()
