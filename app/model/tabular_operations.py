# -IMPORTS FOR INTELLISENSE/TYPE HINTING ONLY-#
from typing import TYPE_CHECKING, Callable
if TYPE_CHECKING:
    from uuid import UUID
    from model.article_managers import Article
    from model.database import TableDBManager
    from pandas import DataFrame

from collections import defaultdict
import logging
import numpy as np
import os
from skimage.measure import label, regionprops
import pandas as pd

from model.file_io import download_supp, extract_dfs
from model.database import processed_df_to_db
from model.article_managers import ProcessedTable

# TODO heavily overdue a readability pass

def cols_to_sorted_lists(table: 'ProcessedTable', df: 'DataFrame'):
    rows = []

    for mapping in table.mappings:
        id_ranges = {}
        for range_start, col in mapping.ids:
            if col not in id_ranges:
                id_ranges[col] = range_start
            else:
                id_ranges[col] = min(id_ranges[col], range_start)
                
        order = "Ascending" if mapping.ascending else "Descending"

        num_ranges = {}
        for range_start, col in mapping.values:
            if col not in num_ranges:
                num_ranges[col] = range_start + 1  # Adjusting for header skip
            else:
                num_ranges[col] = min(num_ranges[col], range_start + 1)

        for num_col, range_start in num_ranges.items():
            column_name = df.columns[num_col]

            temp_df = df.iloc[range_start:].copy()

            try:
                temp_df[column_name] = pd.to_numeric(temp_df[column_name], errors='coerce')
            except Exception as e:
                logging.error(f"Error in cols_to_sorted_lists: {e}")
                pass
            
            sorted_df = temp_df.sort_values(by=column_name, ascending=mapping.ascending)

            for id_col, id_range_start in id_ranges.items():
                table_id = table.id.replace(' ', '_')
                id_col_name = str(df.iat[id_range_start, id_col]).replace(' ', '_')
                value_col_name = str(df.iat[range_start - 1, num_col]).replace(' ', '_')

                list_label = f"[[{table_id}]]__[[{id_col_name}]]__[[{value_col_name}]]__[[{order}]]"
                sorted_list = sorted_df.iloc[:, id_col].tolist()
                list_str = ','.join(map(str, sorted_list))

                rows.append([list_label, list_str])

    result_df = pd.DataFrame(rows, columns=['list_label', 'ranked_list'])
    return result_df


def parse_tables(
        selected_articles: 'list[Article]',
        db_manager: 'TableDBManager',
        should_stop: 'bool',
        callback: 'Callable' = None
):
    tables_metadatas = defaultdict(list)
    for index, article in enumerate(selected_articles):
        processed_table_ids: 'list[tuple[str, UUID]]' = []
        for file in article.supp_files:
            if should_stop:
                break

            if not file.checked:
                continue

            try:
                fname = download_supp(file.url, should_stop)
                if fname is None:
                    continue

                data = extract_dfs(fname, should_stop)

                for sheetname, df in data.items():
                    if should_stop:
                        break

                    if df.empty:
                        continue

                    binary_rep = np.array(df.notnull().astype("int"))
                    labeled = label(binary_rep)
                    region_bboxes = [
                        region.bbox for region in regionprops(labeled)]

                    non_table_bboxes = []  # Step 2

                    for i, box in enumerate(region_bboxes):
                        if should_stop:
                            break

                        minr, minc, maxr, maxc = box
                        if maxr - minr <= 1 or maxc - minc <= 1:
                            non_table_bboxes.append(box)
                            continue

                        other_bboxes = [o_box for j, o_box
                                        in enumerate(region_bboxes)
                                        if i != j]

                        is_contained = any(
                            _is_contained(box, o_box)
                            for o_box in other_bboxes)

                        if is_contained:
                            continue

                        data: 'DataFrame' = df.iloc[minr:maxr, minc:maxc]
                        base_name = os.path.splitext(fname)[0]
                        unique_id = f"{base_name}_{sheetname}_Table{i}"
                        tables_metadatas[unique_id].extend(
                            _find_closest_metadata(box, non_table_bboxes, df)
                        )

                        processed_df_to_db(
                            db_manager,
                            article.pmc_id,
                            article.title,
                            article.url,
                            file.url,
                            data,
                            unique_id,
                            []
                        )

                        processed_table_ids.append((unique_id, file.id))

            except Exception as e:
                logging.exception(
                    f"An error occurred while processing {file.url}: {str(e)}")

        if callback:
            progress = int(100 * (index + 1) / len(selected_articles))
            callback(article, progress, processed_table_ids)


def _find_closest_metadata(
    box: 'tuple',
    non_table_bboxes: 'list[tuple]',
    df: 'DataFrame'
) -> 'list[DataFrame]':
    minr, minc, maxr, maxc = box
    center_x, center_y = (minr + maxr) / 2, (minc + maxc) / 2

    closest_metadata = None
    min_distance = float('inf')

    for nt_box in non_table_bboxes:
        nt_minr, nt_minc, nt_maxr, nt_maxc = nt_box
        nt_center_x, nt_center_y = (
            nt_minr + nt_maxr) / 2, (nt_minc + nt_maxc) / 2

        distance = np.sqrt((center_x - nt_center_x)**2 +
                           (center_y - nt_center_y)**2)

        if distance < min_distance:
            min_distance = distance
            closest_metadata = df.iloc[nt_minr:nt_maxr, nt_minc:nt_maxc]

    return [closest_metadata] if closest_metadata is not None else []


def _is_contained(box: 'tuple', o_box: 'tuple') -> 'bool':
    return (o_box[0] <= box[0] and
            o_box[1] <= box[1] and
            o_box[2] >= box[2] and
            o_box[3] >= box[3])
