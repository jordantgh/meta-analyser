from PyQt5.QtWidgets import QApplication
from threading import Thread, Event
import argparse
import sys
import toml
from datetime import datetime
import os

from model.model import Model
from views.view import View
from controller.controller import Controller

def debug(shutdown_flag, **kwargs):
    def start_repl(local_vars):
        from ptpython.repl import embed

        try:
            while not shutdown_flag.is_set():
                embed(globals(), local_vars)
        except SystemExit:
            pass
        finally:
            sys.exit()

    repl_thread = Thread(target=start_repl, args=(kwargs,))
    repl_thread.start()

def main():
    parser = argparse.ArgumentParser(description="Run the PyQt5 app.")
    parser.add_argument('--debug', action='store_true', help='Activate the debug terminal')
    args = parser.parse_args()
    
    date = datetime.today().strftime('%Y-%m-%d')
    config = toml.load('config.toml')
    
    path = os.path.expanduser(config['data']['path'])
    db_path = f"{path}/db/{date}"
    saves_path = f"{path}/saves/{date}"

    os.makedirs(db_path, exist_ok=True)
    os.makedirs(saves_path, exist_ok=True)

    app = QApplication([])

    model = Model(db_path, saves_path)
    view = View()
    controller = Controller(model, view)

    shutdown_flag = Event()

    if args.debug:
        debug(shutdown_flag, model=model, view=view, controller=controller, app=app)

    view.show()
    app.exec_()

    shutdown_flag.set()

if __name__ == "__main__":
    main()
