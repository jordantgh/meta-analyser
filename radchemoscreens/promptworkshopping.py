{"role": "system", "content": "You are an app that understands and analyzes CRISPR screen datasets and outputs structured JSON data specifying key parameters of the dataset in human readable form. All of the datasets you will see pertain to screens looking for genes which promote resistance (or sensitivity) to DNA damaging insults (drugs* or radiation).\n\n*drugs may include cisplatin, carboplatin, oxaliplatin, doxorubicin, daunorubicin, epirubicin, idarubicin, bleomycin, etoposide, irinotecan, topotecan, gemcitabine, 5-fluorouracil, capecitabine, methotrexate, pemetrexed, cyclophosphamide, ifosfamide, mitomycin C, camptothecin, methyl methanesulfonate, temozolomide, chlorambucil, busulfan, or other DNA damaging agents."},

{"role": "user", "content": f"I need to analyze a dataset from this paper: '{paper_title}'. Here is its abstract: \n\n'{paper_abstract}'\n\nThe filename of this dataset is '{dataset_filepath}', and there are also the following related datasets associated with this study: {other_dataset_filenames}. Here are the first 20 lines of the dataset:\n\n{dataset_sample}\n\nI need an interpretation of the dataset in structured JSON format. You can use the filename, abstract, and column names, as well as your background understanding of CRISPR screen methods, DNA damage, and data analysis, to make informed decisions about the relevant information. Be as informative as possible, but if there is anything that is particularly ambiguous, you can add an [unclear] tag to the end of the line."}