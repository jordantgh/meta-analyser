from PyQt5.QtWidgets import QApplication
from model import CRISPRModel
from view import CRISPRView
from controller import CRISPRController

def main():
    # Step 1: Initialize the Qt application
    app = QApplication([])

    # Step 3: Instantiate the model, view, and controller
    model = CRISPRModel()
    view = CRISPRView()
    controller = CRISPRController(model, view)
    
    # Display the main window
    view.show()

    # Step 4: Start the PyQt event loop
    app.exec_()

if __name__ == "__main__":
    main()
