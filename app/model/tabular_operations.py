import logging
import numpy as np
import os
from skimage.measure import label, regionprops

from model.file_io import download_supp, extract_dfs
from model.database import processed_df_to_db

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

                    for i, bbox1 in enumerate(region_bboxes):
                        if should_stop:
                            break

                        minr1, minc1, maxr1, maxc1 = bbox1
                        if maxr1 - minr1 <= 1 or maxc1 - minc1 <= 1:
                            continue
                        
                        for i, bbox1 in enumerate(region_bboxes):
                            if should_stop:
                                break

                            other_bboxes = [bbox2 for j, bbox2
                                            in enumerate(region_bboxes)
                                            if i != j]
                            
                            is_contained = any(
                                _is_contained(bbox1, bbox2)
                                for bbox2 in other_bboxes)

                            if is_contained:
                                continue
                        
                        region = df.iloc[minr1:maxr1, minc1:maxc1]
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
            callback(article, processed_table_ids, progress)


def _is_contained(bbox1, bbox2):
    return (bbox2[0] <= bbox1[0] and
            bbox2[1] <= bbox1[1] and
            bbox2[2] >= bbox1[2] and
            bbox2[3] >= bbox1[3])