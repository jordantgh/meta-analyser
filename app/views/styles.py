# Using a python file rather than qss for convenience during build
style="""
/* General Styles */
QWidget {
    font-family: "Segoe UI", "Ubuntu", "Helvetica Neue", "Lucida Grande", "Arial", sans-serif;
    font-size: 11pt;
}

QPushButton:focus {
    border: 1px solid blue;
    border-radius: 4px;
}

QLineEdit:focus {
    border: 1px solid blue;
    border-radius: 4px;
}

QListWidget:focus {
    border: 2px solid blue;
}
"""