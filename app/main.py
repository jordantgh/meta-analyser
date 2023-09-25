from PyQt5.QtWidgets import QApplication
from model import Model
from view import View
from controller import Controller
from threading import Thread


def debug(**kwargs):
    """Start an interactive REPL with the given local variables."""
    def start_repl(local_vars):
        from ptpython.repl import embed
        embed(globals(), local_vars)
        
    repl_thread = Thread(target=start_repl, args=(kwargs,))
    repl_thread.start()


def main():
    app = QApplication([])

    model = Model()
    view = View()
    controller = Controller(model, view)

    debug(model=model, view=view, controller=controller, app=app)

    view.show()
    app.exec_()


if __name__ == "__main__":
    main()
