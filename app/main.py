import os
os.environ["XLWINGS_LICENSE_KEY"] = "noncommercial"
os.environ['NLTK_DATA'] = os.path.join(".", "nltk_data")

from PyQt5.QtWidgets import QApplication
import toml
from datetime import datetime
from model.model import Model
from views.view import View
from controller.controller import Controller
import sys

def get_app_data_path():
    if getattr(sys, 'frozen', False):
        # Running as compiled binary
        root = os.path.dirname(sys.executable)
    else:
        # Running as Python script
        root = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(root, 'data')

def main():    
    app = QApplication([])

    # Define default configuration
    default_config = {
        'data': {
            'path': get_app_data_path()
        }
    }

    config_path = "config.toml"
    if not os.path.exists(config_path):
        # Save the default config to a file
        with open(config_path, 'w') as f:
            toml.dump(default_config, f)

    # Load the config file
    config = toml.load(config_path)

    date = datetime.today().strftime('%Y-%m-%d')
    path = os.path.expanduser(config['data']['path'])
    
    db_tpath = os.path.join(path, 'db', 'temp')
    db_ppath = os.path.join(path, 'db', date)
    saves_path = os.path.join(path, 'saves', date)

    os.makedirs(db_tpath, exist_ok=True)
    os.makedirs(db_ppath, exist_ok=True)
    os.makedirs(saves_path, exist_ok=True)

    model = Model(db_tpath, db_ppath, saves_path)
    view = View()
    _ = Controller(model, view)

    view.show()
    app.exec_()


if __name__ == "__main__":
    main()
