# Meta-Analyser

## Overview

Meta-Analyser is a tool designed to automate data search from literature. It performs searches of PubMed Central open access articles (using the normal PMC search syntax) and automatically retrieves supplementary files. It includes a processing workflow that attempts to parse out discrete data frames from often inconsistently formatted .xlsx and .csv files provided with academic publications. processing of scientific articles. Users can preview data files and select the ones they want to process. Once the processing is complete, the user again has the opportunity to further refine the data by selecting the specific tables and features they want to keep. The data is stored in a local SQLite database.

## Features

- Search for scientific articles based on a query
- Download supplementary files linked to the articles
- Process tables from supplementary files and store them in a local SQLite database

## Dependencies

- Python 3.11
- PyQt5
- pandas
- SQLAlchemy
- scikit-image
- requests

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/your-username/meta-analyser.git
    ```

2. Navigate to the project directory:

    ```bash
    cd meta-analyser
    ```

3. Install the required packages using poetry:

    ```bash
    poetry install
    ```

4. Activate the virtual environment:

    ```bash
    poetry shell
    ```

## Usage

1. Run the application:

    ```bash
    python app/main.py
    ```

2. The GUI will appear. Enter your search query into the "Enter Query" field.

3. Click the "Search" button to start searching for articles.

4. Once articles are displayed, select the ones you want to proceed with and click "Proceed".

5. The application will download the supplementary files and process any tables found, storing them in a local SQLite database (`tables.db`).

## Code Structure

- `app/main.py`: Entry point of the application.
- `app/model.py`: Contains the logic for data processing.
- `app/view.py`: Defines the UI components of the application.
- `app/controller.py`: Connects the model and the view.
- `app/search_for_papers.py`: Logic for searching scientific articles (not shown in the code snippet).