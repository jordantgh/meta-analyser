import logging
import numpy as np
import os
from skimage.measure import label, regionprops

from model.file_io import download_supp, extract_dfs
from model.database import processed_df_to_db

# TODO #31 idea: use the contentful regions discarded by the parser to
# add metadata; these are often descriptions/titles of the tables.
# Could then use pythagoras to find the closest table and set the metadata
# on that table.


# TODO heavily overdue a readability pass
def parse_tables(selected_articles, db_manager, should_stop, callback=None):
    for index, article in enumerate(selected_articles):
        processed_table_ids = []
        for file in article.supp_files:
            if should_stop:
                break

            if not file.checked:
                continue

            try:
                fname = download_supp(file.url)
                if fname is None:
                    continue

                data = extract_dfs(fname)

                for sheetname, df in data.items():
                    if should_stop:
                        break

                    if df.empty:
                        continue

                    binary_rep = np.array(df.notnull().astype("int"))
                    labeled = label(binary_rep)
                    region_bboxes = [
                        region.bbox for region in regionprops(labeled)]

                    for i, box in enumerate(region_bboxes):
                        if should_stop:
                            break

                        minr, minc, maxr, maxc = box
                        if maxr - minr <= 1 or maxc - minc <= 1:
                            continue

                        other_bboxes = [o_box for j, o_box
                                        in enumerate(region_bboxes)
                                        if i != j]

                        is_contained = any(
                            _is_contained(box, o_box)
                            for o_box in other_bboxes)

                        if is_contained:
                            continue

                        region = df.iloc[minr:maxr, minc:maxc]
                        base_name = os.path.splitext(fname)[0]
                        unique_id = f"{base_name}_{sheetname}_Table{i}"
                        processed_df_to_db(
                            db_manager, unique_id, file.id, region)

                        processed_table_ids.append((unique_id, file.id))

            except Exception as e:
                logging.exception(
                    f"An error occurred while processing {file.url}: {str(e)}")

        if callback:
            progress = int(100 * (index + 1) / len(selected_articles))
            callback(article, progress, processed_table_ids)


def _is_contained(box, o_box):
    return (o_box[0] <= box[0] and
            o_box[1] <= box[1] and
            o_box[2] >= box[2] and
            o_box[3] >= box[3])
