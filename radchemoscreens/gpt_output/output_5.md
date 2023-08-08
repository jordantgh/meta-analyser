# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1I_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['LORD MCF10A_TP53/RB1mut Olaparib CRISPRn', 'LORD MCF10A_TP53mut Olaparib CRISPRn', 'LORD MCF10A_TP53mut Talazoparib CRISPRn', 'LORD MCF10A_TP53mut Olaparib CRISPRi', 'LORD MCF10A_TP53mut Talazoparib CRISPRi', 'Olivieri RPE1 Olaparib CRISRPn', 'Zimmerman HELA Olaparib CRISPRn', 'DeWeirdt A375 Talazoparib CRISPRn']

condition_variables: ['Study', 'Cell line', 'Drug', 'CRISPR type']

sample:
|    | 0        | 1         | 2           | 3        | 4        | 5             | 6       |
|---:|:---------|:----------|:------------|:---------|:---------|:--------------|:--------|
|  0 | id       | pos_score | pos_p-value | pos_fdr  | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | KEAP1    | 2.87e-10  | 2.59e-07    | 0.002475 | 1        | 3             | 10.566  |
|  2 | C19orf43 | 9.42e-07  | 2.85e-06    | 0.018152 | 2        | 3             | 8.7163  |
|  3 | PIGA     | 2.12e-05  | 7.56e-05    | 0.294554 | 3        | 2             | 4.5552  |
|  4 | BBC3     | 2.24e-05  | 8e-05       | 0.294554 | 4        | 4             | 7.5085  |
|  5 | PRKCSH   | 2.64e-05  | 9.25e-05    | 0.294554 | 5        | 2             | 1.8063  |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_5953_DMSO_Day14_R2_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['DMSO']

condition_variables: ['Drug']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG |  677 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 1248 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   88 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG |  842 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG |  972 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG |  749 |


# ---
dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t2

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin R1', 'Puromycin R2']

condition_variables: ['Drug']

sample:
|    | 0                         | 1            | 2            |
|---:|:--------------------------|:-------------|:-------------|
|  0 | sgRNA                     | Puromycin R1 | Puromycin R2 |
|  1 | A1BG_CATCTTCTTTCACCTGAACG | 1391         | 472          |
|  2 | A1BG_CTCCGGGGAGAACTCCGGCG | 959          | 864          |
|  3 | A1BG_TCTCCATGGTGCATCAGCAC | 13           | 0            |
|  4 | A1BG_TGGAAGTCCACTCCACTCAG | 1937         | 874          |
|  5 | A1CF_ACAGGAAGAATTCAGTTATG | 961          | 740          |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1A_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['Olaparib', 'Talazoparib', 'CRISPRn', 'CRISPRi']

condition_variables: ['Drug', 'CRISPR type']

sample:
|    | 0       | 1                                                 | 2                                         | 3                                         | 4                          | 5                                          | 6                                            | 7                             | 8                                       | 9                                         | 10                                        | 11                         | 12                                       | 13                                         | 14                                           | 15                            |
|---:|:--------|:--------------------------------------------------|:------------------------------------------|:------------------------------------------|:---------------------------|:-------------------------------------------|:---------------------------------------------|:------------------------------|:----------------------------------------|:------------------------------------------|:------------------------------------------|:---------------------------|:-----------------------------------------|:-------------------------------------------|:---------------------------------------------|:------------------------------|
|  0 | GeneID  | CRISPRn Viability Z score (olap and talaz screen) | CRISPRn Drug Effect Z score (olap screen) | CRISPRn Olaparib sgRNA Drug Effect Zcount | CRISPRn Olaparib Gene Rank | CRISPRn Drug Effect Z score (talaz screen) | CRISPRn Talazoparib sgRNA Drug Effect Zcount | CRISPRn Talazoparib Gene Rank | CRISPRi Viability Z score (olap screen) | CRISPRi Drug Effect Z score (olap screen) | CRISPRi Olaparib sgRNA Drug Effect Zcount | CRISPRi Olaparib Gene Rank | CRISPRi Viability Z score (talaz screen) | CRISPRi Drug Effect Z score (talaz screen) | CRISPRi Talazoparib sgRNA Drug Effect Zcount | CRISPRi Talazoparib Gene Rank |
|  1 | EME1    | -1.18066632086832                                 | -5.79777277823934                         | 2                                         | 145.0281356421             | -6.88296766788751                          | 4                                            | 92.189511686447               | -1.04162440274998                       | -5.67188088363195                         | 4                                         | 69.2429995486647           | -1.55784195358798                        | -28.4590476173672                          | 0                                            | 1.25992104989487              |
|  2 | PSMC3IP | -0.978568839197009                                | -11.6888114784934                         | 4                                         | 7.86222418262669           | -13.1931800606307                          | 4                                            | 4.48140474655716              | -5.78049844148705                       | -9.87470478092304                         | 5                                         | 7.73061405251592           | -7.77629543970368                        | -28.3143829575241                          | 3                                            | 1.81712059283214              |
|  3 | MUS81   | -1.56819104215569                                 | -8.3983499013803                          | 3                                         | 27.4470606881676           | -9.17049683624061                          | 3                                            | 20.8561598693187              | -4.31827391718099                       | -9.81232958093168                         | 5                                         | 12.2538513504568           | -5.06825987036802                        | -24.4660198054462                          | 1                                            | 3.30192724889463              |
|  4 | XRCC1   | -0.475994504948211                                | -6.83989924511046                         | 4                                         | 57.1174351203621           | -10.7341472048973                          | 5                                            | 6.83990378670679              | 0.124870862776963                       | -3.9538092680304                          | 1                                         | 270.776475406793           | -0.166748749560336                       | -22.5125138292494                          | 0                                            | 6                             |
|  5 | LIG1    | -1.67118342488582                                 | -13.2487238483912                         | 5                                         | 7.48887238721851           | -11.0493873440622                          | 5                                            | 9.35609523735209              | -1.62330129540551                       | -12.5744032604674                         | 4                                         | 3.30192724889463           | -1.19348469717757                        | -18.7400640806795                          | 0                                            | 6.83990378670679              |


# ---
dataset_filename: caiCooperationATMFanconi2020_32075772_NIHMS1563617-supplement-2_CP0041_20170705_compat_chip_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

condition_variables: ['None']

sample:
|    | 0                    | 1                     | 2                 |
|---:|:---------------------|:----------------------|:------------------|
|  0 | Barcode Sequence     | Annotated Gene Symbol | Annotated Gene ID |
|  1 | AAAAAAAATCCGGACAATGG | SLC25A24              | 29957             |
|  2 | AAAAAAAGGATGGTGATCAA | FASTKD3               | 79072             |
|  3 | AAAAAAATGACATTACTGCA | BCAS2                 | 10286             |
|  4 | AAAAAAATGTCAGTCGAGTG | GPR18                 | 2841              |
|  5 | AAAAAACACAAGCAAGACCG | ZNF470                | 388566            |


# ---
dataset_filename: chenGenomicLandscapeSensitivity2023_37251921_Table_1_ATO_normalized-total.gene_summa_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)']

conditions: ['ATO']

condition_variables: ['Drug']

sample:
|    | 0        | 1         | 2           | 3        | 4        | 5             | 6       |
|---:|:---------|:----------|:------------|:---------|:---------|:--------------|:--------|
|  0 | id       | pos_score | pos_p-value | pos_fdr  | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | KEAP1    | 2.87e-10  | 2.59e-07    | 0.002475 | 1        | 3             | 10.566  |
|  2 | C19orf43 | 9.42e-07  | 2.85e-06    | 0.018152 | 2        | 3             | 8.7163  |
|  3 | PIGA     | 2.12e-05  | 7.56e-05    | 0.294554 | 3        | 2             | 4.5552  |
|  4 | BBC3     | 2.24e-05  | 8e-05       | 0.294554 | 4        | 4             | 7.5085  |
|  5 | PRKCSH   | 2.64e-05  | 9.25e-05    | 0.294554 | 5        | 2             | 1.8063  |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_5953_DMSO_Day14_R2_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['DMSO']

condition_variables: ['Drug']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG |  677 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 1248 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   88 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG |  842 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG |  972 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG |  749 |


# ---
dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t2

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

condition_variables: ['Drug']

sample:
|    | 0                         | 1            | 2            |
|---:|:--------------------------|:-------------|:-------------|
|  0 | sgRNA                     | Puromycin R1 | Puromycin R2 |
|  1 | A1BG_CATCTTCTTTCACCTGAACG | 1391         | 472          |
|  2 | A1BG_CTCCGGGGAGAACTCCGGCG | 959          | 864          |
|  3 | A1BG_TCTCCATGGTGCATCAGCAC | 13           | 0            |
|  4 | A1BG_TGGAAGTCCACTCCACTCAG | 1937         | 874          |
|  5 | A1CF_ACAGGAAGAATTCAGTTATG | 961          | 740          |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1A_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['CRISPRn Olaparib', 'CRISPRn Talazoparib', 'CRISPRi Olaparib', 'CRISPRi Talazoparib']

condition_variables: ['CRISPR type', 'Drug']

sample:
|    | 0       | 1                                                 | 2                                         | 3                                         | 4                          | 5                                          | 6                                            | 7                             | 8                                       | 9                                         | 10                                        | 11                         | 12                                       | 13                                         | 14                                           | 15                            |
|---:|:--------|:--------------------------------------------------|:------------------------------------------|:------------------------------------------|:---------------------------|:-------------------------------------------|:---------------------------------------------|:------------------------------|:----------------------------------------|:------------------------------------------|:------------------------------------------|:---------------------------|:-----------------------------------------|:-------------------------------------------|:---------------------------------------------|:------------------------------|
|  0 | GeneID  | CRISPRn Viability Z score (olap and talaz screen) | CRISPRn Drug Effect Z score (olap screen) | CRISPRn Olaparib sgRNA Drug Effect Zcount | CRISPRn Olaparib Gene Rank | CRISPRn Drug Effect Z score (talaz screen) | CRISPRn Talazoparib sgRNA Drug Effect Zcount | CRISPRn Talazoparib Gene Rank | CRISPRi Viability Z score (olap screen) | CRISPRi Drug Effect Z score (olap screen) | CRISPRi Olaparib sgRNA Drug Effect Zcount | CRISPRi Olaparib Gene Rank | CRISPRi Viability Z score (talaz screen) | CRISPRi Drug Effect Z score (talaz screen) | CRISPRi Talazoparib sgRNA Drug Effect Zcount | CRISPRi Talazoparib Gene Rank |
|  1 | EME1    | -1.18066632086832                                 | -5.79777277823934                         | 2                                         | 145.0281356421             | -6.88296766788751                          | 4                                            | 92.189511686447               | -1.04162440274998                       | -5.67188088363195                         | 4                                         | 69.2429995486647           | -1.55784195358798                        | -28.4590476173672                          | 0                                            | 1.25992104989487              |
|  2 | PSMC3IP | -0.978568839197009                                | -11.6888114784934                         | 4                                         | 7.86222418262669           | -13.1931800606307                          | 4                                            | 4.48140474655716              | -5.78049844148705                       | -9.87470478092304                         | 5                                         | 7.73061405251592           | -7.77629543970368                        | -28.3143829575241                          | 3                                            | 1.81712059283214              |
|  3 | MUS81   | -1.56819104215569                                 | -8.3983499013803                          | 3                                         | 27.4470606881676           | -9.17049683624061                          | 3                                            | 20.8561598693187              | -4.31827391718099                       | -9.81232958093168                         | 5                                         | 12.2538513504568           | -5.06825987036802                        | -24.4660198054462                          | 1                                            | 3.30192724889463              |
|  4 | XRCC1   | -0.475994504948211                                | -6.83989924511046                         | 4                                         | 57.1174351203621           | -10.7341472048973                          | 5                                            | 6.83990378670679              | 0.124870862776963                       | -3.9538092680304                          | 1                                         | 270.776475406793           | -0.166748749560336                       | -22.5125138292494                          | 0                                            | 6                             |
|  5 | LIG1    | -1.67118342488582                                 | -13.2487238483912                         | 5                                         | 7.48887238721851           | -11.0493873440622                          | 5                                            | 9.35609523735209              | -1.62330129540551                       | -12.5744032604674                         | 4                                         | 3.30192724889463           | -1.19348469717757                        | -18.7400640806795                          | 0                                            | 6.83990378670679              |


# ---
dataset_filename: caiCooperationATMFanconi2020_32075772_NIHMS1563617-supplement-2_CP0041_20170705_compat_chip_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Lookup/Reference/Metadata

statistic_aggregation: NA

metrics: ['NA']

conditions: ['NA']

condition_variables: ['NA']

sample:
|    | 0                    | 1                     | 2                 |
|---:|:---------------------|:----------------------|:------------------|
|  0 | Barcode Sequence     | Annotated Gene Symbol | Annotated Gene ID |
|  1 | AAAAAAAATCCGGACAATGG | SLC25A24              | 29957             |
|  2 | AAAAAAAGGATGGTGATCAA | FASTKD3               | 79072             |
|  3 | AAAAAAATGACATTACTGCA | BCAS2                 | 10286             |
|  4 | AAAAAAATGTCAGTCGAGTG | GPR18                 | 2841              |
|  5 | AAAAAACACAAGCAAGACCG | ZNF470                | 388566            |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1I_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['LORD_MCF10ATP53/RB1mut_Olaparib_CRISPRn', 'LORD_MCF10ATP53mut_Olaparib_CRISPRn', 'LORD_MCF10A TP53mut_Talazoparib_CRISPRn', 'LORD_MCF10ATP53mut_Olaparib_CRISPRi', 'LORD_MCF10ATP53mut_Talazoparib_CRISPRi', 'Olivieri_RPE1_Olaparib_CRISRPn', 'Zimmerman_HELA_Olaparib_CRISPRn', 'DeWeirdt_A375_Talazoparib_CRISPRn']

condition_variables: ['Study', 'Cell line', 'Drug', 'CRISPR type']

sample:
|    | 0       | 1                                       | 2                                   | 3                                       | 4                                   | 5                                      | 6                              | 7                               | 8                                 |
|---:|:--------|:----------------------------------------|:------------------------------------|:----------------------------------------|:------------------------------------|:---------------------------------------|:-------------------------------|:--------------------------------|:----------------------------------|
|  0 | GeneID  | LORD_MCF10ATP53/RB1mut_Olaparib_CRISPRn | LORD_MCF10ATP53mut_Olaparib_CRISPRn | LORD_MCF10A TP53mut_Talazoparib_CRISPRn | LORD_MCF10ATP53mut_Olaparib_CRISPRi | LORD_MCF10ATP53mut_Talazoparib_CRISPRi | Olivieri_RPE1_Olaparib_CRISRPn | Zimmerman_HELA_Olaparib_CRISPRn | DeWeirdt_A375_Talazoparib_CRISPRn |
|  1 | A1BG    | 0.0072091                               | 0.3717573                           | 0.4623897                               | -0.5694143                          | 0.1544734                              | 1.2778839                      | -1.1010888                      | -0.5029758                        |
|  2 | A1CF    | -2.7335572                              | -0.9436333                          | 2.0481784                               | 0.1533515                           | 0.6418993                              | -1.0936161                     | 0.6558407                       | -0.3588095                        |
|  3 | A2M     | 0.0738308                               | 0.5541413                           | 0.0683911                               | -0.2250607                          | 1.3635369                              | -0.8166494                     | 0.8628098                       | 0.6884674                         |
|  4 | A2ML1   | 1.5033446                               | 1.1912513                           | 0.94766                                 | 0.1072929                           | 1.7982553                              | -1.0152515                     | 0.8943272                       | -0.0798228                        |
|  5 | A3GALT2 | 0.6887889                               | 0.1351586                           | -0.1016254                              | 0.166886                            | 0.7562855                              | 0.2436562                      | nan                             | 1.2573963                         |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_6MP_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)', 'MAGeCK statistics (negative selection)']

conditions: ['6MP']

condition_variables: ['Drug']

sample:
|    | 0       | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8          | 9           | 10      | 11       | 12            | 13      |
|---:|:--------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:-----------|:------------|:--------|:---------|:--------------|:--------|
|  0 | id      | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score  | pos_p-value | pos_fdr | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | HPRT1   | 6   | 0.99989   | 0.99991     | 1       | 20371    | 0             | 2.5776  | 3.2452e-21 | 2.4189e-07  | 0.00019 | 1        | 6             | 2.5776  |
|  2 | SLC43A3 | 6   | 0.99996   | 0.99997     | 1       | 20395    | 0             | 1.6773  | 2.0736e-17 | 2.4189e-07  | 0.00019 | 2        | 6             | 1.6773  |
|  3 | NUDT5   | 6   | 1         | 1           | 1       | 20466    | 0             | 1.7439  | 3.4419e-17 | 2.4189e-07  | 0.00019 | 3        | 6             | 1.7439  |
|  4 | CSTF3   | 6   | 1         | 1           | 1       | 20465    | 0             | 0.99532 | 1.1683e-12 | 2.4189e-07  | 0.00019 | 4        | 6             | 0.99532 |
|  5 | EIF2B3  | 6   | 1         | 1           | 1       | 20464    | 0             | 1.0834  | 1.3033e-12 | 2.4189e-07  | 0.00019 | 5        | 6             | 1.0834  |


# ---
dataset_filename: goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 5_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol [MOD]

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Cisplatin', 'DMSO']

condition_variables: ['Drug']

sample:
|    | 0                                                                       | 1          | 2          | 3          | 4          | 5     | 6     | 7     | 8     |
|---:|:------------------------------------------------------------------------|:-----------|:-----------|:-----------|:-----------|:------|:------|:------|:------|
|  0 | Supplementary Table 5. Raw counts from the cisplatin resistance screen. | nan        | nan        | nan        | nan        | nan   | nan   | nan   | nan   |
|  1 | sgRNA                                                                   | Cisplatin1 | Cisplatin2 | Cisplatin3 | Cisplatin4 | DMSO1 | DMSO2 | DMSO3 | DMSO4 |
|  2 | EHMT1-HGLibA_14685                                                      | 123        | 228        | 147        | 16         | 224   | 234   | 208   | 89    |
|  3 | RALGPS2-HGLibA_40172                                                    | 5          | 2          | 7          | 167        | 15    | 155   | 36    | 90    |
|  4 | HIST2H2BE-HGLibA_21534                                                  | 9          | 492        | 1          | 49         | 28    | 0     | 72    | 91    |
|  5 | SGCA-HGLibA_43795                                                       | 173        | 87         | 273        | 82         | 142   | 111   | 174   | 156   |


# ---
dataset_filename: liuSpindleAssemblyCheckpoint2019_30862715_NIHMS1524158-supplement-5_Cisplatin_DMF_gene_summary.txt_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (negative selection)']

conditions: ['Cisplatin']

condition_variables: ['Drug']

sample:
|    | 0      | 1                        | 2                          | 3                       |
|---:|:-------|:-------------------------|:---------------------------|:------------------------|
|  0 | id     | negative selection score | negative selection p-value | negative selection rank |
|  1 | GLIS1  | 2.6779e-05               | 0.00014369                 | 1                       |
|  2 | CNTRL  | 2.9913e-05               | 0.00015625                 | 2                       |
|  3 | EFNB2  | 4.0944e-05               | 0.0002065                  | 3                       |
|  4 | CCDC54 | 6.045e-05                | 0.00027821                 | 4                       |
|  5 | RGAG4  | 8.0334e-05               | 0.00036929                 | 5                       |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1J_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['CRISPRn Olaparib']

condition_variables: ['CRISPR type', 'Drug']

sample:
|    | 0       | 1                         | 2                               | 3                              |
|---:|:--------|:--------------------------|:--------------------------------|:-------------------------------|
|  0 | GeneID  | CRISPRn Viability Z-score | CRISPRn Olaparib Effect Z-score | CRISPRn Olaparib sgRNA Z-count |
|  1 | A1BG    | -0.670949539095226        | 0.0944227963684928              | 0                              |
|  2 | A1CF    | 0.758921722245498         | -4.03972344728084               | 1                              |
|  3 | A2M     | -0.21371060796212         | 0.169043715504459               | 0                              |
|  4 | A2ML1   | 0.529894455033704         | 1.83532457638626                | 0                              |
|  5 | A3GALT2 | -1.09668301945057         | 0.895692189469412               | 0                              |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Vincristine_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['Vincristine']

condition_variables: ['Drug']

sample:
|    | 0       | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8         | 9           | 10      | 11       | 12            | 13      |
|---:|:--------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:----------|:------------|:--------|:---------|:--------------|:--------|
|  0 | id      | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score | pos_p-value | pos_fdr | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | TP53    | 6   | 0.97305   | 0.97497     | 1       | 19574    | 1             | 0.93454 | 2.44e-12  | 2.42e-07    | 0.00045 | 1        | 5             | 0.93454 |
|  2 | TRMT112 | 6   | 1         | 1           | 1       | 20458    | 0             | 0.9559  | 3.77e-11  | 2.42e-07    | 0.00045 | 2        | 5             | 0.9559  |
|  3 | RRP9    | 6   | 1         | 1           | 1       | 20466    | 0             | 0.67563 | 4.62e-11  | 2.42e-07    | 0.00045 | 3        | 6             | 0.67563 |
|  4 | CRCP    | 6   | 0.99904   | 0.99904     | 1       | 20295    | 0             | 0.91053 | 1.72e-10  | 2.42e-07    | 0.00045 | 4        | 5             | 0.91053 |
|  5 | NOP16   | 6   | 1         | 1           | 1       | 20463    | 0             | 0.99505 | 1.82e-10  | 2.42e-07    | 0.00045 | 5        | 6             | 0.99505 |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Daunorubicin_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['Daunorubicin']

condition_variables: ['Drug']

sample:
|    | 0     | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8          | 9           | 10       | 11       | 12            | 13      |
|---:|:------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:-----------|:------------|:---------|:---------|:--------------|:--------|
|  0 | id    | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score  | pos_p-value | pos_fdr  | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | HINFP | 6   | 1         | 1           | 1       | 20466    | 0             | 1.8969  | 5.9965e-16 | 2.4189e-07  | 0.000248 | 1        | 6             | 1.8969  |
|  2 | TOP2B | 6   | 0.95021   | 0.95021     | 1       | 18965    | 0             | 2.4823  | 4.504e-13  | 2.4189e-07  | 0.000248 | 2        | 5             | 2.4823  |
|  3 | CSTF3 | 6   | 1         | 1           | 1       | 20465    | 0             | 1.1544  | 3.9642e-12 | 2.4189e-07  | 0.000248 | 3        | 6             | 1.1544  |
|  4 | CRCP  | 6   | 1         | 1           | 1       | 20463    | 0             | 1.6789  | 2.6661e-11 | 2.4189e-07  | 0.000248 | 4        | 6             | 1.6789  |
|  5 | IARS2 | 6   | 1         | 1           | 1       | 20464    | 0             | 1.5182  | 4.2423e-11 | 2.4189e-07  | 0.000248 | 5        | 6             | 1.5182  |


# ---
dataset_filename: biaynaLossAbasicSite2021_33788831_pbio.3001176.s014_lfc_all_time_points_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['fold-change/log fold-change']

conditions: ['A549 TP53-/- time 9', 'A549 TP53-/- time 12', 'A549 TP53-/- time 15', 'A549 TP53 WT time 9', 'A549 TP53 WT time 12', 'A549 TP53 WT time 15', 'LXF289 DOX-IC25 time 5', 'LXF289 DOX-IC25 time 10', 'LXF289 DOX-IC25 time 15', 'LXF289 DOX-IC50 time 5', 'LXF289 DOX-IC50 time 10', 'LXF289 DOX-IC50 time 15']

condition_variables: ['Cell line', 'Time', 'Drug']

sample:
|    | 0     | 1            | 2         | 3         | 4            | 5         | 6         | 7               | 8         | 9         | 10              | 11        | 12        |
|---:|:------|:-------------|:----------|:----------|:-------------|:----------|:----------|:----------------|:----------|:----------|:----------------|:----------|:----------|
|  0 | nan   | A549 TP53-/- | nan       | nan       | A549 TP53 WT | nan       | nan       | LXF289 DOX-IC25 | nan       | nan       | LXF289 DOX-IC50 | nan       | nan       |
|  1 | Gene  | time 9       | time 12   | time 15   | time 9       | time 12   | time 15   | time 5          | time 10   | time 15   | time 5          | time 10   | time 15   |
|  2 | A1BG  | 0.24838      | -0.26098  | 0.099242  | 0.35132      | -0.025552 | -0.01677  | -0.07502        | -0.073809 | -0.22727  | 0.051521        | 0.041967  | -0.069046 |
|  3 | A1CF  | -0.006276    | 0.0038166 | -0.002651 | -0.25883     | 0.07213   | 0.060406  | -0.0064161      | -0.084975 | 0.057137  | -0.1384         | -0.14075  | 0.1142    |
|  4 | A2M   | 0.12738      | -0.1451   | -0.12855  | 0.19655      | 0.072992  | -0.012698 | 0.084314        | 0.28355   | -0.013266 | -0.15669        | 0.0010223 | -0.51692  |
|  5 | A2ML1 | -0.0030069   | -0.039999 | -0.064889 | 0.13849      | 0.11629   | 0.07048   | -0.1096         | -0.022392 | 0.2551    | -0.19769        | 0.051254  | 0.24108   |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_6398_DMSO_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['DMSO']

condition_variables: ['Drug']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG |  474 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 1175 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   82 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG |  573 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG | 1212 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG | 1374 |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_6377_Library_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

condition_variables: ['None']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG | 1021 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 1485 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   26 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG | 1665 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG | 1068 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG |  802 |


# ---
dataset_filename: ramakerPooledCRISPRScreening2021_34049503_12885_2021_8388_MOESM1_ESM_TableS3A_t1

header_row: 3

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['P1 Gemcitabine', 'P1 Irinotecan', 'P1 Oxaliplatin', 'P1 5FU', 'B3 Gemcitabine', 'B3 Irinotecan', 'B3 Oxaliplatin', 'B3 5FU']

condition_variables: ['Cell line', 'Drug']

sample:
|    | 0                                                                                                                         | 1    | 2                 | 3           | 4            | 5                       | 6                | 7                  | 8                      | 9               | 10                | 11                      | 12               | 13                 | 14              | 15          | 16          | 17                      | 18               | 19                 | 20                     | 21              | 22                | 23                      | 24               | 25                 | 26              | 27          | 28          |
|---:|:--------------------------------------------------------------------------------------------------------------------------|:-----|:------------------|:------------|:-------------|:------------------------|:-----------------|:-------------------|:-----------------------|:----------------|:------------------|:------------------------|:-----------------|:-------------------|:----------------|:------------|:------------|:------------------------|:-----------------|:-------------------|:-----------------------|:----------------|:------------------|:------------------------|:-----------------|:-------------------|:----------------|:------------|:------------|
|  0 | Supplemental Table S3A. Activation (SAM) Screen results. L2FC Sum for individual drug-cell line combinations and combined | nan  | nan               | nan         | nan          | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         |
|  1 | P1= Panc1,B3=BxPC3                                                                                                        | nan  | nan               | nan         | nan          | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         |
|  2 | For each gene in each screen we calculated a L2FC sum, a p-value and FDR                                                  | nan  | nan               | nan         | nan          | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan         | nan         |
|  3 | Gene_ID                                                                                                                   | nan  | Combined_L2FC_Sum | Combined_P  | Combined_FDR | P1_Gemcitabine_L2FC_Sum | P1_Gemcitabine_P | P1_Gemcitabine_FDR | P1_Irinotecan_L2FC_Sum | P1_Irinotecan_P | P1_Irinotecan_FDR | P1_Oxaliplatin_L2FC_Sum | P1_Oxaliplatin_P | P1_Oxaliplatin_FDR | P1_5FU_L2FC_Sum | P1_5FU_P    | P1_5FU_FDR  | B3_Gemcitabine_L2FC_Sum | B3_Gemcitabine_P | B3_Gemcitabine_FDR | B3_Irinotecan_L2FC_Sum | B3_Irinotecan_P | B3_Irinotecan_FDR | B3_Oxaliplatin_L2FC_Sum | B3_Oxaliplatin_P | B3_Oxaliplatin_FDR | B3_5FU_L2FC_Sum | B3_5FU_P    | B3_5FU_FDR  |
|  4 | NM_000014                                                                                                                 | A2M  | -0.18554717       | 0.130953394 | 0.350599116  | 1.30817191              | 0.02168621       | 0.129813642        | -0.093118004           | 0.050648391     | 0.16322417        | 1.565240464             | 0.004759342      | 0.048339743        | 0.31002847      | 0.171653606 | 0.394650874 | -0.087805172            | 0.594800367      | 0.890740617        | -0.170216254           | 0.62710939      | 0.916172815       | -0.147466101            | 0.669739761      | 0.911610932        | -0.336701153    | 0.744820208 | 0.948290624 |
|  5 | NM_000015                                                                                                                 | NAT2 | -0.615084931      | 0.226631071 | 0.380142049  | 0.404445504             | 0.218215912      | 0.461666725        | -2.04439538            | 0.231691291     | 0.428193874       | -1.161697916            | 0.827891999      | 0.944957403        | 0.341308066     | 0.160581291 | 0.379965492 | 0.163152039             | 0.252576833      | 0.693515925        | 0.306888492            | 0.131873646     | 0.603273276       | 0.520695716             | 0.021178454      | 0.285921004        | -0.246447184    | 0.665934009 | 0.921367956 |


# ---
dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts RP-6306 resistance_t1

header_row: 1

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['TKOv3_FT282.CCNE1oe_C3_NT_T6A', 'TKOv3_FT282.CCNE1oe_C3_NT_T6B', 'TKOv3_FT282.CCNE1oe_C4_NT_T6A', 'TKOv3_FT282.CCNE1oe_C4_NT_T6B', 'TKOv3_FT282.CCNE1oe_C3_RP6306_T21A', 'TKOv3_FT282.CCNE1oe_C3_RP6306_T21B', 'TKOv3_FT282.CCNE1oe_C4_RP6306_T21A', 'TKOv3_FT282.CCNE1oe_C4_RP6306_T21B']

condition_variables: ['Cell line', 'Drug']

sample:
|    | 0                                                                                             | 1        | 2                             | 3                             | 4                             | 5                             | 6                                  | 7                                  | 8                                  | 9                                  |
|---:|:----------------------------------------------------------------------------------------------|:---------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:-----------------------------------|:-----------------------------------|:-----------------------------------|:-----------------------------------|
|  0 | Read counts for FT282-hTERT p53-R175H CCNE1-high RP-6306 Resistance Screen for Clones 3 and 4 | nan      | nan                           | nan                           | nan                           | nan                           | nan                                | nan                                | nan                                | nan                                |
|  1 | sgRNA                                                                                         | Gene     | TKOv3_FT282.CCNE1oe_C3_NT_T6A | TKOv3_FT282.CCNE1oe_C3_NT_T6B | TKOv3_FT282.CCNE1oe_C4_NT_T6A | TKOv3_FT282.CCNE1oe_C4_NT_T6B | TKOv3_FT282.CCNE1oe_C3_RP6306_T21A | TKOv3_FT282.CCNE1oe_C3_RP6306_T21B | TKOv3_FT282.CCNE1oe_C4_RP6306_T21A | TKOv3_FT282.CCNE1oe_C4_RP6306_T21B |
|  2 | chr11:134201957-134201976_GLB1L2_+                                                            | GLB1L2   | 1737                          | 1521                          | 1549                          | 1406                          | 708                                | 1100                               | 1178                               | 843                                |
|  3 | chr6:26022026-26022045_HIST1H4A_-                                                             | HIST1H4A | 2052                          | 1379                          | 2108                          | 1631                          | 1735                               | 975                                | 1723                               | 853                                |
|  4 | chr12:120436384-120436403_CCDC64_+                                                            | CCDC64   | 425                           | 391                           | 496                           | 476                           | 146                                | 506                                | 337                                | 92                                 |
|  5 | chr8:73979678-73979697_SBSPON_-                                                               | SBSPON   | 894                           | 800                           | 1048                          | 780                           | 956                                | 1193                               | 1100                               | 1100                               |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_6397_Puromycin_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

condition_variables: ['Drug']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG |  373 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG |  814 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   10 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG |  344 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG |  880 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG | 1049 |


# ---
dataset_filename: dengIdentifyingCDC7Synergistic2023_36725843_41420_2023_1315_MOESM8_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)', 'MAGeCK statistics (negative selection)']

conditions: ['d7_chemo', 'd7_control']

condition_variables: ['Drug']

sample:
|    | 0        | 1                         | 2                                 | 3                             | 4                              | 5                                   | 6                                 | 7                             | 8                              | 9                                   | 10                                   |
|---:|:---------|:--------------------------|:----------------------------------|:------------------------------|:-------------------------------|:------------------------------------|:----------------------------------|:------------------------------|:-------------------------------|:------------------------------------|:-------------------------------------|
|  0 | id       | d7_chemo__d7_control__num | d7_chemo__d7_control__neg_p-value | d7_chemo__d7_control__neg_fdr | d7_chemo__d7_control__neg_rank | d7_chemo__d7_control__neg_goodsgrna | d7_chemo__d7_control__pos_p-value | d7_chemo__d7_control__pos_fdr | d7_chemo__d7_control__pos_rank | d7_chemo__d7_control__pos_goodsgrna | d7_chemo__d7_control__log2FoldChange |
|  1 | TRMT12   | 3                         | 0.29732                           | 0.996322                      | 6050                           | 1                                   | 5.3276e-06                        | 0.113861                      | 1                              | 2                                   | 9.1497                               |
|  2 | C11orf16 | 1                         | 0.99868                           | 0.999829                      | 21344                          | 0                                   | 0.0013627                         | 0.956958                      | 68                             | 1                                   | 9.0236                               |
|  3 | FAM208B  | 1                         | 0.99693                           | 0.999632                      | 21312                          | 0                                   | 0.0030884                         | 0.956958                      | 152                            | 1                                   | 8.8196                               |
|  4 | APOBR    | 1                         | 0.996                             | 0.999471                      | 21297                          | 0                                   | 0.0039862                         | 0.956958                      | 196                            | 1                                   | 8.7494                               |
|  5 | WFDC13   | 1                         | 0.99597                           | 0.999471                      | 21296                          | 0                                   | 0.004027                          | 0.956958                      | 197                            | 1                                   | 8.7472                               |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1C_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Z-score']

conditions: ['CRISPRi Olaparib', 'CRISPRi Talazoparib']

condition_variables: ['CRISPR type', 'Drug']

sample:
|    | 0      | 1               | 2                     | 3                                       | 4                                         | 5                                        | 6                                          | 7             |
|---:|:-------|:----------------|:----------------------|:----------------------------------------|:------------------------------------------|:-----------------------------------------|:-------------------------------------------|:--------------|
|  0 | GeneID | Ensembl gene ID | sgRNA                 | CRISPRi Viability Z score (olap screen) | CRISPRi Drug Effect Z score (olap screen) | CRISPRi Viability Z score (talaz screen) | CRISPRi Drug Effect Z score (talaz screen) | Essentiality  |
|  1 | A1BG   | ENSG00000121410 | A1BG_+_58864367.23-P2 | 0.125489761                             | 0.470384573                               | 0.317818471                              | -0.551194082                               | non-essential |
|  2 | A1BG   | ENSG00000121410 | A1BG_-_58864840.23-P2 | 0.956319208                             | -0.077996751                              | 1.340699837                              | -0.413110529                               | non-essential |
|  3 | A1BG   | ENSG00000121410 | A1BG_-_58864822.23-P2 | -0.451517775                            | 0.449761143                               | -0.524998017                             | 0.637647537                                | non-essential |
|  4 | A1BG   | ENSG00000121410 | A1BG_+_58858549.23-P1 | 0.805701784                             | -0.589471123                              | 0.717304988                              | 0.011958931                                | non-essential |
|  5 | A1BG   | ENSG00000121410 | A1BG_-_58864360.23-P2 | -1.056138343                            | 0.070177808                               | -1.200402162                             | -0.628214011                               | non-essential |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_L-asparaginase_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['L-asparaginase']

condition_variables: ['Drug']

sample:
|    | 0      | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8         | 9           | 10      | 11       | 12            | 13      |
|---:|:-------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:----------|:------------|:--------|:---------|:--------------|:--------|
|  0 | id     | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score | pos_p-value | pos_fdr | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | DHODH  | 6   | 0.98227   | 0.98211     | 1       | 19328    | 0             | 2.9941  | 1.09e-11  | 2.42e-07    | 0.00016 | 1        | 5             | 2.9941  |
|  2 | G6PD   | 6   | 1         | 1           | 1       | 20466    | 0             | 1.7795  | 1.37e-11  | 2.42e-07    | 0.00016 | 2        | 6             | 1.7795  |
|  3 | CAD    | 6   | 1         | 1           | 1       | 20465    | 0             | 1.8716  | 1.42e-11  | 2.42e-07    | 0.00016 | 3        | 6             | 1.8716  |
|  4 | PMPCA  | 6   | 1         | 1           | 1       | 20464    | 0             | 1.4795  | 2.08e-11  | 2.42e-07    | 0.00016 | 4        | 6             | 1.4795  |
|  5 | IMPDH2 | 6   | 1         | 1           | 1       | 20463    | 0             | 1.6361  | 6.2e-11   | 2.42e-07    | 0.00016 | 5        | 6             | 1.6361  |


# ---
dataset_filename: masoudiGenomescaleCRISPRCas92019_31844142_41598_2019_55893_MOESM2_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)']

conditions: ['Gemcitabine']

condition_variables: ['Drug']

sample:
|    | 0            | 1          | 2        | 3         | 4       | 5            |
|---:|:-------------|:-----------|:---------|:----------|:--------|:-------------|
|  0 | Gene         | # Hairpins | NES      | Gene rank | p-value | p-value rank |
|  1 | SH3D21       | 6          | 0.003022 | 1         | 0.0001  | 2            |
|  2 | DNAJB12      | 6          | 0.003166 | 2         | 0.0001  | 3            |
|  3 | hsa-mir-4443 | 4          | 0.003691 | 3         | 0.0001  | 1            |
|  4 | COA3         | 6          | 0.005052 | 4         | 0.0003  | 5            |
|  5 | ATAD2        | 6          | 0.007434 | 5         | 0.0004  | 6            |


# ---
dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Raw count']

conditions: ['Eto R1', 'Eto R2', 'DMSO']

condition_variables: ['Drug']

sample:
|    | 0     | 1      | 2      | 3    |
|---:|:------|:-------|:-------|:-----|
|  0 | genes | Eto R1 | Eto R2 | DMSO |
|  1 | A1BG  | 2032   | 286    | 677  |
|  2 | A1BG  | 934    | 497    | 1248 |
|  3 | A1BG  | 0      | 8      | 88   |
|  4 | A1BG  | 2544   | 726    | 842  |
|  5 | A1CF  | 1191   | 236    | 972  |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_6396_Puromycin_Day0_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

condition_variables: ['Drug']

sample:
|    | 0                         |    1 |
|---:|:--------------------------|-----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG |  608 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 1154 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   31 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG | 1008 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG |  852 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG | 1040 |


# ---
dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t3

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Library R1', 'Library R2']

condition_variables: ['Library']

sample:
|    | 0                         | 1          | 2          |
|---:|:--------------------------|:-----------|:-----------|
|  0 | sgRNA                     | Library R1 | Library R2 |
|  1 | A1BG_CATCTTCTTTCACCTGAACG | 1373       | 1004       |
|  2 | A1BG_CTCCGGGGAGAACTCCGGCG | 1188       | 1496       |
|  3 | A1BG_TCTCCATGGTGCATCAGCAC | 21         | 14         |
|  4 | A1BG_TGGAAGTCCACTCCACTCAG | 2090       | 1652       |
|  5 | A1CF_ACAGGAAGAATTCAGTTATG | 1101       | 1053       |


# ---
dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts TKOv2 CCNE1 SL_t1

header_row: 1

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#2 T0', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#2 T15_repA', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#2 T15_repB', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#2 T18_repA', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#2 T18_repB', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#21 T0', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#21 T15_repA', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#21 T15_repB', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#21 T18_repA', 'TKOv2 RPE1.hTERT.Cas9.dTP53 CCNE1oe Clone#21 T18_repB', 'TKOv2 RPE1.hTERT.Cas9.dTP53_T0', 'TKOv2 RPE1.hTERT.Cas9.dTP53_WT_T18_repA', 'TKOv2 RPE1.hTERT.Cas9.dTP53_T18_repB']

condition_variables: ['Cell line', 'Time point', 'Replicate']

sample:
|    | 0                                   | 1     | 2                                              | 3                                                    | 4                                                    | 5                                                    | 6                                                    | 7                                               | 8                                                     | 9                                                     | 10                                                    | 11                                                    | 12                             | 13                                      | 14                                   |
|---:|:------------------------------------|:------|:-----------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:-------------------------------|:----------------------------------------|:-------------------------------------|
|  0 | Read Counts for screens using TKOv2 | nan   | nan                                            | nan                                                  | nan                                                  | nan                                                  | nan                                                  | nan                                             | nan                                                   | nan                                                   | nan                                                   | nan                                                   | nan                            | nan                                     | nan                                  |
|  1 | sgRNA                               | Gene  | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T15_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T15_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T15_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T15_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T18_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_WT_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_T18_repB |
|  2 | chr1:100111914-100111933_PALMD_-    | PALMD | 516                                            | 242                                                  | 166                                                  | 325                                                  | 283                                                  | 1106                                            | 220                                                   | 521                                                   | 305                                                   | 117                                                   | 1340                           | 690                                     | 676                                  |
|  3 | chr1:100133222-100133241_PALMD_+    | PALMD | 101                                            | 35                                                   | 31                                                   | 52                                                   | 37                                                   | 129                                             | 72                                                    | 99                                                    | 50                                                    | 35                                                    | 167                            | 117                                     | 62                                   |
|  4 | chr1:100152265-100152284_PALMD_+    | PALMD | 808                                            | 290                                                  | 134                                                  | 370                                                  | 255                                                  | 703                                             | 194                                                   | 902                                                   | 218                                                   | 179                                                   | 1316                           | 705                                     | 586                                  |
|  5 | chr1:100154717-100154736_PALMD_+    | PALMD | 425                                            | 239                                                  | 133                                                  | 170                                                  | 201                                                  | 476                                             | 122                                                   | 453                                                   | 145                                                   | 153                                                   | 439                            | 421                                     | 277                                  |


# ---
dataset_filename: funkeGenomescaleCRISPRScreen2023_37024667_41416_2023_2247_MOESM3_ESM_7_NGS_Cand_list_t2

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: RefSeq ID

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

condition_variables: ['None']

sample:
|    |   0 | 1         | 2            | 3                     | 4           |
|---:|----:|:----------|:-------------|:----------------------|:------------|
|  0 | nan | Gene Name | Refseq ID    | Target sequence sgRNA | Read counts |
|  1 | nan | CYP2R1    | NM_024514    | AGAAGGTAGTTGTCCCGAAG  | 4248699     |
|  2 | nan | LEKR1     | NM_001193283 | GGAGCGTTGACGACAGACGT  | 4084952     |
|  3 | nan | CCL1      | NM_002981    | GAGTTTACCACATGGTCGCT  | 2230750     |
|  4 | nan | LHB       | NM_000894    | CATGTGCACCTCTCGCCCCC  | 1223512     |
|  5 | nan | WBP2NL    | NM_152613    | AGCGCGAAGGGGTGAGTGAT  | 1066313     |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Methotrexate_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)', 'MAGeCK statistics (negative selection)']

conditions: ['Methotrexate']

condition_variables: ['Drug']

sample:
|    | 0       | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8          | 9           | 10      | 11       | 12            | 13      |
|---:|:--------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:-----------|:------------|:--------|:---------|:--------------|:--------|
|  0 | id      | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score  | pos_p-value | pos_fdr | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | SLC19A1 | 6   | 0.97481   | 0.97475     | 1       | 19354    | 0             | 4.8868  | 6.637e-21  | 2.4189e-07  | 8e-05   | 1        | 5             | 4.8868  |
|  2 | NAA60   | 6   | 0.99972   | 0.99973     | 1       | 20257    | 0             | 4.4439  | 1.0139e-17 | 2.4189e-07  | 8e-05   | 2        | 6             | 4.4439  |
|  3 | PPP2CA  | 5   | 0.99927   | 0.99929     | 1       | 20200    | 0             | 3.2935  | 3.6789e-09 | 2.4189e-07  | 8e-05   | 45       | 5             | 3.2935  |
|  4 | CNOT4   | 6   | 1         | 1           | 1       | 20405    | 0             | 3.0542  | 3.6357e-16 | 2.4189e-07  | 8e-05   | 3        | 6             | 3.0542  |
|  5 | HNRNPM  | 6   | 0.98809   | 0.98801     | 1       | 19724    | 0             | 2.9752  | 8.8629e-13 | 2.4189e-07  | 8e-05   | 8        | 5             | 2.9752  |


# ---
dataset_filename: ramakerPooledCRISPRScreening2021_34049503_12885_2021_8388_MOESM1_ESM_TableS3B_t1

header_row: 3

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score']

conditions: ['P1 Gemcitabine', 'P1 Oxaliplatin', 'P1 Irinotecan', 'P1 5FU', 'B3 Gemcitabine', 'B3 Oxaliplatin', 'B3 Irinotecan', 'B3 5FU']

condition_variables: ['Cell line', 'Drug']

sample:
|    | 0                                                                                                            | 1                  | 2                  | 3                 | 4                       | 5                  | 6                  | 7                       | 8                  | 9                  | 10                     | 11                 | 12                | 13                | 14                 | 15                | 16                      | 17                | 18                 | 19                      | 20                | 21                 | 22                     | 23                | 24                | 25                | 26                 | 27                |
|---:|:-------------------------------------------------------------------------------------------------------------|:-------------------|:-------------------|:------------------|:------------------------|:-------------------|:-------------------|:------------------------|:-------------------|:-------------------|:-----------------------|:-------------------|:------------------|:------------------|:-------------------|:------------------|:------------------------|:------------------|:-------------------|:------------------------|:------------------|:-------------------|:-----------------------|:------------------|:------------------|:------------------|:-------------------|:------------------|
|  0 | Supplemental Table S3B. Knockout (Gecko) Screen results. L2FC Sum for individual drug-cell line combinations | nan                | nan                | nan               | nan                     | nan                | nan                | nan                     | nan                | nan                | nan                    | nan                | nan               | nan               | nan                | nan               | nan                     | nan               | nan                | nan                     | nan               | nan                | nan                    | nan               | nan               | nan               | nan                | nan               |
|  1 | P1 = Panc1: B3=BxPC3                                                                                         | nan                | nan                | nan               | nan                     | nan                | nan                | nan                     | nan                | nan                | nan                    | nan                | nan               | nan               | nan                | nan               | nan                     | nan               | nan                | nan                     | nan               | nan                | nan                    | nan               | nan               | nan               | nan                | nan               |
|  2 | For each gene in each screen we calculated a L2FC sum, a p-value and FDR                                     | nan                | nan                | nan               | nan                     | nan                | nan                | nan                     | nan                | nan                | nan                    | nan                | nan               | nan               | nan                | nan               | nan                     | nan               | nan                | nan                     | nan               | nan                | nan                    | nan               | nan               | nan               | nan                | nan               |
|  3 | nan                                                                                                          | Combined_L2FC_Sum  | Combined_P         | Combined_FDR      | P1_Gemcitabine_L2FC_Sum | P1_Gemcitabine_P   | P1_Gemcitabine_FDR | P1_Oxaliplatin_L2FC_Sum | P1_Oxaliplatin_P   | P1_Oxaliplatin_FDR | P1_Irinotecan_L2FC_Sum | P1_Irinotecan_P    | P1_Irinotecan_FDR | P1_5FU_L2FC_Sum   | P1_5FU_P           | P1_5FU_FDR        | B3_Gemcitabine_L2FC_Sum | B3_Gemcitabine_P  | B3_Gemcitabine_FDR | B3_Oxaliplatin_L2FC_Sum | B3_Oxaliplatin_P  | B3_Oxaliplatin_FDR | B3_Irinotecan_L2FC_Sum | B3_Irinotecan_P   | B3_Irinotecan_FDR | B3_5FU_L2FC_Sum   | B3_5FU_P           | B3_5FU_FDR        |
|  4 | A1BG                                                                                                         | -0.904586216999054 | 0.314775208036283  | 0.701825580880729 | -0.904586216999054      | 0.705751098984387  | 0.999965452086164  | -0.392715133913581      | 0.72360982484397   | 0.970318452715167  | -1.4862571781002       | 0.28926700246282   | 0.999898342971155 | -2.04481219892761 | 0.618499289440127  | 0.99991789536423  | -0.378696386660878      | 0.747680849435469 | 0.993359561111788  | -0.00544005991251131    | 0.394241943897826 | 0.962196700272364  | -0.607283237222274     | 0.781932864584173 | 0.999628192829801 | -1.11520453166738 | 0.950561172993587  | 0.999687840103799 |
|  5 | A1CF                                                                                                         | 0.351737973373145  | 0.0193697286062979 | 0.650760987670613 | 0.351737973373145       | 0.0927310393613259 | 0.922223015712754  | 0.517615043804666       | 0.0444151399234951 | 0.618721214655904  | -0.107825095499609     | 0.0502095521796057 | 0.807361064319172 | 0.37677409411936  | 0.0518009297082954 | 0.689169112092335 | -0.684541032563738      | 0.903139463221877 | 0.997791930332779  | -0.362510690205391      | 0.673465770503558 | 0.981022096637038  | -0.378322136192466     | 0.641481478495057 | 0.999628192829801 | 0.510627260795835 | 0.0663233750561576 | 0.885783329272473 |


# ---
dataset_filename: awahGenomeScaleCRISPR2022_35990011_6399_Temozolomide_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Temozolomide']

condition_variables: ['Drug']

sample:
|    | 0                         |   1 |
|---:|:--------------------------|----:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG | 311 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG | 355 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC |   4 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG | 487 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG | 512 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG | 589 |


# ---
dataset_filename: funkeGenomescaleCRISPRScreen2023_37024667_41416_2023_2247_MOESM3_ESM_7_NGS_Cand_list_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: RefSeq ID

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

condition_variables: ['None']

sample:
|    | 0         | 1         | 2                     | 3           |   4 |
|---:|:----------|:----------|:----------------------|:------------|----:|
|  0 | Gene Name | Refseq ID | Target sequence sgRNA | Read counts | nan |
|  1 | TRAP1     | NM_016292 | AGGGCGACGGGCCTTGCGCG  | 8044032     | nan |
|  2 | TK1       | NM_003258 | CGTGCGTCCCTCTGTTTATA  | 7803883     | nan |
|  3 | DDB1      | NM_001923 | GCCTTCGGATGTGGCGGATG  | 6841409     | nan |
|  4 | ENSA      | NM_207047 | TGCTTTGGCGCTGGTTAGTT  | 1590607     | nan |
|  5 | IFNW1     | NM_002177 | CCCAGGTACTCAGGGGGCTG  | 656648      | nan |


# ---
dataset_filename: nguyenGenomewideCRISPRCas92022_36111891_BLOODA_ADV-2022-007934-mmc1_Supplemental Table S0_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['Etoposide', 'Vehicle control']

condition_variables: ['Drug', 'Time']

sample:
|    | 0                                                                                                                                                                                                                          | 1                  | 2             | 3              | 4            | 5                  | 6              | 7               | 8             | 9                  | 10             | 11              | 12            |
|---:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------|:--------------|:---------------|:-------------|:-------------------|:---------------|:----------------|:--------------|:-------------------|:---------------|:----------------|:--------------|
|  0 | Supplemental Table S0. Genome-wide Etoposide score (etoposide vs. vehicle controls) and essential score (vehicle control vs. day 0 replicates) at day 4 (early), day 12 (intermediate), and day 18 (late) using MAGeCK-RRA | nan                | nan           | nan            | nan          | nan                | nan            | nan             | nan           | nan                | nan            | nan             | nan           |
|  1 | Gene                                                                                                                                                                                                                       | Day4.Etop.Score    | Day4.Etop.FDR | Day4.Ess.Score | Day4.Ess.FDR | Day12.Etop.Score   | Day12.Etop.FDR | Day12.Ess.Score | Day12.Ess.FDR | Day18.Etop.Score   | Day18.Etop.FDR | Day18.Ess.Score | Day18.Ess.FDR |
|  2 | A1BG                                                                                                                                                                                                                       | -0.982757915452354 | 0.972306      | 0.20568        | 0.853545     | 1.7519285990808    | 0.752928       | -1.7533         | 0.099903      | -1.33446591156667  | 0.933154       | -0.036103       | 1             |
|  3 | A1CF                                                                                                                                                                                                                       | 0.341254264638847  | 1             | 0.092881       | 1            | -0.483398528473466 | 1              | 0.08054         | 1             | 1.19622888552981   | 0.960285       | 0.22909         | 1             |
|  4 | A2M                                                                                                                                                                                                                        | -0.621456752838356 | 0.998435      | 0.1778         | 0.970935     | -0.715817505363098 | 1              | 0.56573         | 0.82638       | 1.07435277487943   | 0.979498       | 0.18195         | 0.97963       |
|  5 | A2ML1                                                                                                                                                                                                                      | 0.726373793727041  | 1             | -0.10209       | 0.999976     | 0.3131939780091    | 1              | 0.16833         | 1             | -0.995205889611288 | 0.978192       | 0.13114         | 1             |


# ---
dataset_filename: goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 6_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['fold-change/log fold-change']

conditions: ['Cisplatin']

condition_variables: ['Drug']

sample:
|    | 0                                                                | 1           | 2           | 3           | 4           | 5            | 6           | 7           | 8           | 9          | 10   | 11         | 12         | 13                | 14                    |
|---:|:-----------------------------------------------------------------|:------------|:------------|:------------|:------------|:-------------|:------------|:------------|:------------|:-----------|:-----|:-----------|:-----------|:------------------|:----------------------|
|  0 | Supplementary Table 6. Normalized screen data and DeSEQ2 output. | nan         | nan         | nan         | nan         | nan          | nan         | nan         | nan         | nan        | nan  | nan        | nan        | nan               | nan                   |
|  1 | Gene                                                             | baseMean1   | Log2_FC1    | pval1       | baseMean2   | Log2_FC2     | pval2       | baseMean3   | Log2_FC3    | pval3      | chi  | combined_p | combined_q | Meidan_Log2_FC    | Combined_mean_Log2_FC |
|  2 | MSH2                                                             | 957.2566405 | 5.67843972  | 1.83e-11    | 1975.317814 | 5.116799119  | 1.77e-11    | 1072.351948 | 4.822651573 | 2.38e-07   | 187  | 1.1e-37    | 1.83e-33   | 5.11679911917568  | 5.25058214863253      |
|  3 | MLH1                                                             | 1825.532851 | 3.686081025 | 7.95e-08    | 581.8766207 | 3.253182804  | 0.000102254 | 67.24306563 | 2.459707092 | 0.01856056 | 85.2 | 3e-16      | 4.97e-12   | 3.25318280400865  | 3.21759138327125      |
|  4 | FAM89B                                                           | 34.07272074 | 0.309519502 | 1           | 35.94595749 | -4.138685456 | 2.42e-06    | 418.9378725 | 3.69554062  | 5.49e-07   | 78.9 | 6.03e-15   | 9.98e-11   | 0.309519502289633 | 2.2481295828933       |
|  5 | XPC                                                              | 242.0108639 | 3.590859153 | 0.000204781 | 510.0933263 | 2.508409082  | 0.004336734 | 3303.070108 | 3.92785658  | 1.8e-06    | 78.4 | 7.64e-15   | 1.27e-10   | 3.59085915261611  | 3.45762272302659      |


# ---
dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1B_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Z-score']

conditions: ['CRISPRn Viability (olap and talaz screen)', 'CRISPRn Drug Effect (olap screen)', 'CRISPRn Drug Effect (talaz screen)']

condition_variables: ['CRISPR type', 'Drug']

sample:
|    | 0      | 1               | 2                                               | 3                                                 | 4                                         | 5                                          | 6             |
|---:|:-------|:----------------|:------------------------------------------------|:--------------------------------------------------|:------------------------------------------|:-------------------------------------------|:--------------|
|  0 | GeneID | Ensembl gene ID | sgRNA                                           | CRISPRn Viability Z score (olap and talaz screen) | CRISPRn Drug Effect Z score (olap screen) | CRISPRn Drug Effect Z score (talaz screen) | Essentiality  |
|  1 | A1BG   | ENSG00000121410 | A1BG_CCDS12976.1_ex3_19:58862927-58862950:-_5-1 | 0.668991                                          | 0.705699                                  | -0.505318                                  | non-essential |
|  2 | A1BG   | ENSG00000121410 | A1BG_CCDS12976.1_ex5_19:58864367-58864390:-_5-5 | -1.158484                                         | 0.467419                                  | 0.022966                                   | non-essential |
|  3 | A1BG   | ENSG00000121410 | A1BG_CCDS12976.1_ex4_19:58863655-58863678:+_5-2 | -0.090289                                         | 0.260061                                  | 0.814947                                   | non-essential |
|  4 | A1BG   | ENSG00000121410 | A1BG_CCDS12976.1_ex4_19:58863697-58863720:-_5-3 | -1.123212                                         | 0.215305                                  | 0.754094                                   | non-essential |
|  5 | A1BG   | ENSG00000121410 | A1BG_CCDS12976.1_ex4_19:58863866-58863889:+_5-4 | -2.072579                                         | -0.373944                                 | 0.446271                                   | non-essential |


# ---
dataset_filename: xuGenomewideCRISPRScreen2019_31792210_41467_2019_13420_MOESM4_ESM_reads out counts_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Survive', 'Input']

condition_variables: ['Condition']

sample:
|    | 0      | 1         | 2         | 3       | 4       |
|---:|:-------|:----------|:----------|:--------|:--------|
|  0 | GeneID | Survive_1 | Survive_2 | Input_1 | Input_2 |
|  1 | A1BG   | 0         | 0         | 0       | 0       |
|  2 | A1BG   | 396       | 219       | 153     | 154     |
|  3 | A1BG   | 286       | 362       | 226     | 262     |
|  4 | A1BG   | 325       | 237       | 275     | 266     |
|  5 | A1CF   | 128       | 116       | 129     | 132     |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Maphosphamide_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['Maphosphamide']

condition_variables: ['Drug']

sample:
|    | 0      | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8          | 9           | 10      | 11       | 12            | 13      |
|---:|:-------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:-----------|:------------|:--------|:---------|:--------------|:--------|
|  0 | id     | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score  | pos_p-value | pos_fdr | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | TP53   | 6   | 1         | 1           | 1       | 20466    | 0             | 2.1621  | 1.7926e-19 | 2.4189e-07  | 0.00033 | 1        | 6             | 2.1621  |
|  2 | HEATR3 | 6   | 1         | 1           | 1       | 20460    | 0             | 1.2388  | 3.7987e-12 | 2.4189e-07  | 0.00033 | 2        | 6             | 1.2388  |
|  3 | KAT5   | 6   | 0.57038   | 0.71191     | 1       | 14108    | 1             | 1.2079  | 5.6488e-12 | 2.4189e-07  | 0.00033 | 3        | 5             | 1.2079  |
|  4 | HNRNPM | 6   | 0.90892   | 0.91167     | 1       | 18329    | 1             | 1.0584  | 8.3701e-11 | 2.4189e-07  | 0.00033 | 4        | 5             | 1.0584  |
|  5 | PMAIP1 | 6   | 1         | 1           | 1       | 20465    | 0             | 0.84946 | 8.5806e-11 | 2.4189e-07  | 0.00033 | 5        | 6             | 0.84946 |


# ---
dataset_filename: masoudiGenomescaleCRISPRCas92019_31844142_41598_2019_55893_MOESM3_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive selection)']

conditions: ['Gemcitabine']

condition_variables: ['Drug']

sample:
|    | 0        | 1          | 2         | 3         | 4       | 5            |
|---:|:---------|:-----------|:----------|:----------|:--------|:-------------|
|  0 | Gene     | # Hairpins | NES       | Gene rank | p-value | p-value rank |
|  1 | KIAA0196 | 6          | 0.0005764 | 1         | 0.0001  | 4            |
|  2 | SLFN14   | 6          | 0.004255  | 2         | 0.0001  | 2            |
|  3 | MTRNR2L2 | 6          | 0.004413  | 3         | 0.0001  | 1            |
|  4 | DGCR6    | 6          | 0.005611  | 4         | 0.0004  | 6            |
|  5 | TMSB4Y   | 4          | 0.005614  | 5         | 0.0001  | 3            |


# ---
dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts_TKOv3 CCNE1 SL_t1

header_row: 1

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0', 'TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18A', 'TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18B', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T0', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18A', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18B']

condition_variables: ['Cell line', 'Time point']

sample:
|    | 0                                   | 1    | 2                                              | 3                                                | 4                                                | 5                                 | 6                                   | 7                                   |
|---:|:------------------------------------|:-----|:-----------------------------------------------|:-------------------------------------------------|:-------------------------------------------------|:----------------------------------|:------------------------------------|:------------------------------------|
|  0 | Read counts for screens using TKOv3 | nan  | nan                                            | nan                                              | nan                                              | nan                               | nan                                 | nan                                 |
|  1 | sgRNA                               | Gene | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0 | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18A | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18B | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T0 | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18A | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18B |
|  2 | TKOv3_A1BG_exon1_4                  | A1BG | 1501                                           | 1029                                             | 240                                              | 451                               | 366                                 | 440                                 |
|  3 | TKOv3_A1BG_exon3_3                  | A1BG | 1404                                           | 703                                              | 463                                              | 404                               | 352                                 | 404                                 |
|  4 | TKOv3_A1BG_exon4_2                  | A1BG | 535                                            | 401                                              | 266                                              | 146                               | 66                                  | 179                                 |
|  5 | TKOv3_A1BG_exon5_1                  | A1BG | 708                                            | 400                                              | 219                                              | 236                               | 111                                 | 252                                 |


# ---
dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_AraC_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics (positive/negative selection)']

conditions: ['AraC']

condition_variables: ['Drug']

sample:
|    | 0       | 1   | 2         | 3           | 4       | 5        | 6             | 7       | 8          | 9           | 10       | 11       | 12            | 13      |
|---:|:--------|:----|:----------|:------------|:--------|:---------|:--------------|:--------|:-----------|:------------|:---------|:---------|:--------------|:--------|
|  0 | id      | num | neg_score | neg_p-value | neg_fdr | neg_rank | neg_goodsgrna | neg_lfc | pos_score  | pos_p-value | pos_fdr  | pos_rank | pos_goodsgrna | pos_lfc |
|  1 | DCK     | 6   | 1         | 1           | 1       | 20466    | 0             | 3.3873  | 1.1002e-22 | 2.4189e-07  | 0.000177 | 1        | 6             | 3.3873  |
|  2 | SLC29A1 | 6   | 1         | 1           | 1       | 20462    | 0             | 3.2527  | 4.7006e-15 | 2.4189e-07  | 0.000177 | 2        | 6             | 3.2527  |
|  3 | NOP16   | 6   | 1         | 1           | 1       | 20465    | 0             | 1.5237  | 8.6231e-13 | 2.4189e-07  | 0.000177 | 3        | 6             | 1.5237  |
|  4 | IARS2   | 6   | 1         | 1           | 1       | 20464    | 0             | 1.7807  | 2.1605e-12 | 2.4189e-07  | 0.000177 | 4        | 6             | 1.7807  |
|  5 | FARSA   | 6   | 1         | 1           | 1       | 20463    | 0             | 1.2936  | 1.5707e-11 | 2.4189e-07  | 0.000177 | 5        | 6             | 1.2936  |


