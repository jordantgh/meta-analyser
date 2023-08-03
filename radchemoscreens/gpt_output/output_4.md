# ---
article_title: A Whole-genome CRISPR Screen Identifies a Role of MSH2 in Cisplatin-mediated Cell Death in Muscle-invasive Bladder Cancer

dataset_filename: goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 5_t1.csv

header_row: 1

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol [MOD]

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Cisplatin', 'DMSO']

sample:
|    | id       |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:---------|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | KEAP1    |  2.87e-10   |    2.59e-07   |  0.002475 |          1 |               3 |   10.566  |
|  1 | C19orf43 |  9.42e-07   |    2.85e-06   |  0.018152 |          2 |               3 |    8.7163 |
|  2 | PIGA     |  2.12e-05   |    7.56e-05   |  0.294554 |          3 |               2 |    4.5552 |
|  3 | BBC3     |  2.24e-05   |    8e-05      |  0.294554 |          4 |               4 |    7.5085 |
|  4 | PRKCSH   |  2.64e-05   |    9.25e-05   |  0.294554 |          5 |               2 |    1.8063 |
|  5 | NTN1     |  3.62e-05   |    0.00012926 |  0.3529   |          6 |               4 |    6.1809 |
|  6 | UBALD1   |  8.37e-05   |    0.00032095 |  0.721163 |          7 |               3 |    7.3436 |
|  7 | G6PC     |  9.29e-05   |    0.00035255 |  0.721163 |          8 |               3 |    7.3963 |
|  8 | EIF2AK1  |  0.00011096 |    0.00042042 |  0.721163 |          9 |               4 |    7.2147 |
|  9 | OGT      |  0.00011148 |    0.00042249 |  0.721163 |         10 |               3 |    6.1563 |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_5953_DMSO_Day14_R2_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['DMSO']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   677 |
|---:|:----------------------------|------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |  1248 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |    88 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |   842 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |   972 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |   749 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |   509 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |   145 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |   166 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |   130 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |   326 |


# ---
article_title: Ribosomal Protein S11 Influences Glioma Response to TOP2 Poisons

dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t2

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

sample:
|    | sgRNA                      |   Puromycin R1 |   Puromycin R2 |
|---:|:---------------------------|---------------:|---------------:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG  |           1391 |            472 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG  |            959 |            864 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC  |             13 |              0 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG  |           1937 |            874 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG  |            961 |            740 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG  |            793 |            608 |
|  6 | A1CF_CTTCATTTCCCAGCCACCAA  |           1378 |            496 |
|  7 | A1CF_GAATTCAACAATATCAAACC  |            880 |            265 |
|  8 | A2ML1_AAATACTTACTGGTATCGAG |            632 |            224 |
|  9 | A2ML1_ACCCTGGTTACTGATAACAA |            340 |            170 |


# ---
article_title: {{MND1}} and {{PSMC3IP}} Control {{PARP}} Inhibitor Sensitivity in Mitotic Cells

dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1A_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score/Z-ratios']

conditions: ['olaparib', 'talazoparib']

sample:
|    | GeneID   |   CRISPRn Viability Z score (olap and talaz screen) |   CRISPRn Drug Effect Z score (olap screen) |   CRISPRn Olaparib sgRNA Drug Effect Zcount |   CRISPRn Olaparib Gene Rank |   CRISPRn Drug Effect Z score (talaz screen) |   CRISPRn Talazoparib sgRNA Drug Effect Zcount |   CRISPRn Talazoparib Gene Rank |   CRISPRi Viability Z score (olap screen) |   CRISPRi Drug Effect Z score (olap screen) |   CRISPRi Olaparib sgRNA Drug Effect Zcount |   CRISPRi Olaparib Gene Rank |   CRISPRi Viability Z score (talaz screen) |   CRISPRi Drug Effect Z score (talaz screen) |   CRISPRi Talazoparib sgRNA Drug Effect Zcount |   CRISPRi Talazoparib Gene Rank |
|---:|:---------|----------------------------------------------------:|--------------------------------------------:|--------------------------------------------:|-----------------------------:|---------------------------------------------:|-----------------------------------------------:|--------------------------------:|------------------------------------------:|--------------------------------------------:|--------------------------------------------:|-----------------------------:|-------------------------------------------:|---------------------------------------------:|-----------------------------------------------:|--------------------------------:|
|  0 | EME1     |                                          -1.18067   |                                    -5.79777 |                                           2 |                    145.028   |                                     -6.88297 |                                              4 |                        92.1895  |                               -1.04162    |                                    -5.67188 |                                           4 |                     69.243   |                                 -1.55784   |                                     -28.459  |                                              0 |                         1.25992 |
|  1 | PSMC3IP  |                                          -0.978569  |                                   -11.6888  |                                           4 |                      7.86222 |                                    -13.1932  |                                              4 |                         4.4814  |                               -5.7805     |                                    -9.8747  |                                           5 |                      7.73061 |                                 -7.7763    |                                     -28.3144 |                                              3 |                         1.81712 |
|  2 | MUS81    |                                          -1.56819   |                                    -8.39835 |                                           3 |                     27.4471  |                                     -9.1705  |                                              3 |                        20.8562  |                               -4.31827    |                                    -9.81233 |                                           5 |                     12.2539  |                                 -5.06826   |                                     -24.466  |                                              1 |                         3.30193 |
|  3 | XRCC1    |                                          -0.475995  |                                    -6.8399  |                                           4 |                     57.1174  |                                    -10.7341  |                                              5 |                         6.8399  |                                0.124871   |                                    -3.95381 |                                           1 |                    270.776   |                                 -0.166749  |                                     -22.5125 |                                              0 |                         6       |
|  4 | LIG1     |                                          -1.67118   |                                   -13.2487  |                                           5 |                      7.48887 |                                    -11.0494  |                                              5 |                         9.3561  |                               -1.6233     |                                   -12.5744  |                                           4 |                      3.30193 |                                 -1.19348   |                                     -18.7401 |                                              0 |                         6.8399  |
|  5 | MND1     |                                           0.0672105 |                                   -10.8408  |                                           4 |                      8.89592 |                                    -11.9469  |                                              4 |                         9.16566 |                               -0.679034   |                                    -6.0834  |                                           4 |                     30.8252  |                                 -0.689851  |                                     -18.6826 |                                              0 |                         3.91487 |
|  6 | ATM      |                                          -0.405374  |                                    -9.08109 |                                           5 |                     14.518   |                                     -9.48875 |                                              4 |                        16.8122  |                               -0.140151   |                                    -7.87685 |                                           4 |                     33.2379  |                                 -1.66117   |                                     -18.5984 |                                              5 |                        11.8547  |
|  7 | DDX11    |                                          -4.10654   |                                    -4.69825 |                                           2 |                    614.753   |                                     -4.112   |                                              1 |                       523.59    |                               -3.01377    |                                    -5.84625 |                                           3 |                    101.985   |                                 -3.14137   |                                     -18.2982 |                                              1 |                        15.8103  |
|  8 | POLB     |                                           0.615428  |                                    -3.43831 |                                           1 |                    820.249   |                                    -11.8818  |                                              4 |                        11.9722  |                               -0.00176596 |                                    -2.8684  |                                           2 |                    637.384   |                                  0.0289104 |                                     -17.8175 |                                              0 |                         8.96281 |
|  9 | TRAIP    |                                          -3.77465   |                                    -2.76533 |                                           1 |                   1413.61    |                                     -2.77421 |                                              2 |                      1234.45    |                               -8.39735    |                                   -12.4171  |                                           5 |                      5.19249 |                                -10.6089    |                                     -17.4829 |                                              0 |                        13.3887  |


# ---
article_title: Cooperation of the ATM and Fanconi Anemia/BRCA Pathways in Double-Strand Break End Resection

dataset_filename: caiCooperationATMFanconi2020_32075772_NIHMS1563617-supplement-2_CP0041_20170705_compat_chip_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: []

sample:
|    | Barcode Sequence     | Annotated Gene Symbol   |   Annotated Gene ID |
|---:|:---------------------|:------------------------|--------------------:|
|  0 | AAAAAAAATCCGGACAATGG | SLC25A24                |               29957 |
|  1 | AAAAAAAGGATGGTGATCAA | FASTKD3                 |               79072 |
|  2 | AAAAAAATGACATTACTGCA | BCAS2                   |               10286 |
|  3 | AAAAAAATGTCAGTCGAGTG | GPR18                   |                2841 |
|  4 | AAAAAACACAAGCAAGACCG | ZNF470                  |              388566 |
|  5 | AAAAAACAGATGCCACCTGT | CENPC                   |                1060 |
|  6 | AAAAAACCAAACTTTGAAGT | CIAPIN1                 |               57019 |
|  7 | AAAAAACCCGTAGATAGCCT | NDUFA12                 |               55967 |
|  8 | AAAAAACCTGGGCAAAACAG | QTRT2                   |               79691 |
|  9 | AAAAAACCTTCAAAGTACAA | STAP1                   |               26228 |


# ---
article_title: {{MND1}} and {{PSMC3IP}} Control {{PARP}} Inhibitor Sensitivity in Mitotic Cells

dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1I_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Other/ambiguous']

conditions: ['Olaparib', 'Talazoparib']

sample:
|    | GeneID   |   LORD_MCF10ATP53/RB1mut_Olaparib_CRISPRn |   LORD_MCF10ATP53mut_Olaparib_CRISPRn |   LORD_MCF10A TP53mut_Talazoparib_CRISPRn |   LORD_MCF10ATP53mut_Olaparib_CRISPRi |   LORD_MCF10ATP53mut_Talazoparib_CRISPRi |   Olivieri_RPE1_Olaparib_CRISRPn |   Zimmerman_HELA_Olaparib_CRISPRn |   DeWeirdt_A375_Talazoparib_CRISPRn |
|---:|:---------|------------------------------------------:|--------------------------------------:|------------------------------------------:|--------------------------------------:|-----------------------------------------:|---------------------------------:|----------------------------------:|------------------------------------:|
|  0 | A1BG     |                                 0.0072091 |                              0.371757 |                                 0.46239   |                            -0.569414  |                                0.154473  |                        1.27788   |                        -1.10109   |                          -0.502976  |
|  1 | A1CF     |                                -2.73356   |                             -0.943633 |                                 2.04818   |                             0.153352  |                                0.641899  |                       -1.09362   |                         0.655841  |                          -0.35881   |
|  2 | A2M      |                                 0.0738308 |                              0.554141 |                                 0.0683911 |                            -0.225061  |                                1.36354   |                       -0.816649  |                         0.86281   |                           0.688467  |
|  3 | A2ML1    |                                 1.50334   |                              1.19125  |                                 0.94766   |                             0.107293  |                                1.79826   |                       -1.01525   |                         0.894327  |                          -0.0798228 |
|  4 | A3GALT2  |                                 0.688789  |                              0.135159 |                                -0.101625  |                             0.166886  |                                0.756286  |                        0.243656  |                       nan         |                           1.2574    |
|  5 | A4GALT   |                                 0.351977  |                             -0.320726 |                                -0.727516  |                             0.364203  |                                0.449817  |                        0.0032903 |                        -0.48303   |                          -0.252524  |
|  6 | A4GNT    |                                 0.14734   |                              1.05048  |                                 0.967364  |                             1.98439   |                               -0.0201671 |                       -0.969091  |                        -3.53388   |                          -0.955747  |
|  7 | AAAS     |                                 0.175251  |                              0.506268 |                                 0.0908106 |                             1.25947   |                                0.564094  |                        0.653248  |                         0.0103649 |                           1.47713   |
|  8 | AACS     |                                 0.878579  |                              0.937439 |                                -0.140813  |                            -0.0700872 |                                1.89266   |                        0.0333263 |                        -0.273311  |                          -0.339773  |
|  9 | AADAC    |                                -1.42566   |                             -1.13581  |                                 1.31016   |                            -0.0625646 |                                0.594484  |                       -1.16398   |                         0.339138  |                           0.680517  |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_6MP_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['6MP']

sample:
|    | id      |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:--------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | HPRT1   |     6 |     0.99989 |       0.99991 |         1 |      20371 |               0 |   2.5776  |  3.2452e-21 |    2.4189e-07 |   0.00019 |          1 |               6 |   2.5776  |
|  1 | SLC43A3 |     6 |     0.99996 |       0.99997 |         1 |      20395 |               0 |   1.6773  |  2.0736e-17 |    2.4189e-07 |   0.00019 |          2 |               6 |   1.6773  |
|  2 | NUDT5   |     6 |     1       |       1       |         1 |      20466 |               0 |   1.7439  |  3.4419e-17 |    2.4189e-07 |   0.00019 |          3 |               6 |   1.7439  |
|  3 | CSTF3   |     6 |     1       |       1       |         1 |      20465 |               0 |   0.99532 |  1.1683e-12 |    2.4189e-07 |   0.00019 |          4 |               6 |   0.99532 |
|  4 | EIF2B3  |     6 |     1       |       1       |         1 |      20464 |               0 |   1.0834  |  1.3033e-12 |    2.4189e-07 |   0.00019 |          5 |               6 |   1.0834  |
|  5 | GMPPB   |     6 |     0.99964 |       0.99964 |         1 |      20324 |               0 |   1.3734  |  8.241e-12  |    2.4189e-07 |   0.00019 |          6 |               6 |   1.3734  |
|  6 | ATAD3A  |     6 |     0.91954 |       0.91961 |         1 |      18324 |               0 |   1.0824  |  1.0243e-11 |    2.4189e-07 |   0.00019 |          7 |               5 |   1.0824  |
|  7 | ROMO1   |     6 |     1       |       1       |         1 |      20463 |               0 |   1.3279  |  2.6608e-11 |    2.4189e-07 |   0.00019 |          8 |               6 |   1.3279  |
|  8 | EIF2B2  |     6 |     0.99995 |       0.99996 |         1 |      20392 |               0 |   1.446   |  7.0301e-11 |    2.4189e-07 |   0.00019 |          9 |               6 |   1.446   |
|  9 | MIOS    |     6 |     1       |       1       |         1 |      20428 |               0 |   1.0697  |  2.0997e-10 |    2.4189e-07 |   0.00019 |         10 |               6 |   1.0697  |


# ---
article_title: A Whole-genome CRISPR Screen Identifies a Role of MSH2 in Cisplatin-mediated Cell Death in Muscle-invasive Bladder Cancer

dataset_filename: goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 5_t1

header_row: 1

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol [MOD]

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Cisplatin', 'DMSO']

sample:
|    | Supplementary Table 5. Raw counts from the cisplatin resistance screen.    | Unnamed: 1   | Unnamed: 2   | Unnamed: 3   | Unnamed: 4   | Unnamed: 5   | Unnamed: 6   | Unnamed: 7   | Unnamed: 8   |
|---:|:---------------------------------------------------------------------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|
|  0 | sgRNA                                                                      | Cisplatin1   | Cisplatin2   | Cisplatin3   | Cisplatin4   | DMSO1        | DMSO2        | DMSO3        | DMSO4        |
|  1 | EHMT1-HGLibA_14685                                                         | 123          | 228          | 147          | 16           | 224          | 234          | 208          | 89           |
|  2 | RALGPS2-HGLibA_40172                                                       | 5            | 2            | 7            | 167          | 15           | 155          | 36           | 90           |
|  3 | HIST2H2BE-HGLibA_21534                                                     | 9            | 492          | 1            | 49           | 28           | 0            | 72           | 91           |
|  4 | SGCA-HGLibA_43795                                                          | 173          | 87           | 273          | 82           | 142          | 111          | 174          | 156          |
|  5 | OPN1MW-HGLibA_33425                                                        | 42           | 71           | 83           | 111          | 190          | 281          | 197          | 198          |
|  6 | GPR83-HGLibA_20093                                                         | 117          | 53           | 0            | 123          | 76           | 125          | 88           | 5            |
|  7 | PRAMEF2-HGLibA_38306                                                       | 74           | 63           | 3            | 16           | 54           | 78           | 102          | 27           |
|  8 | SDCBP2-HGLibA_43111                                                        | 9            | 11           | 0            | 1            | 14           | 5            | 33           | 63           |
|  9 | COPS7A-HGLibA_10729                                                        | 290          | 575          | 215          | 495          | 719          | 420          | 679          | 553          |


# ---
article_title: Spindle Assembly Checkpoint Inhibition Can Resensitize P53-Null Stem Cells to Cancer Chemotherapy

dataset_filename: liuSpindleAssemblyCheckpoint2019_30862715_NIHMS1524158-supplement-5_Cisplatin_DMF_gene_summary.txt_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Cisplatin']

sample:
|    | id       |   negative selection score |   negative selection p-value |   negative selection rank |
|---:|:---------|---------------------------:|-----------------------------:|--------------------------:|
|  0 | GLIS1    |                 2.6779e-05 |                   0.00014369 |                         1 |
|  1 | CNTRL    |                 2.9913e-05 |                   0.00015625 |                         2 |
|  2 | EFNB2    |                 4.0944e-05 |                   0.0002065  |                         3 |
|  3 | CCDC54   |                 6.045e-05  |                   0.00027821 |                         4 |
|  4 | RGAG4    |                 8.0334e-05 |                   0.00036929 |                         5 |
|  5 | PGLYRP4  |                 0.0001245  |                   0.0005771  |                         6 |
|  6 | PSMC3    |                 0.00013389 |                   0.0006195  |                         7 |
|  7 | ZNF207   |                 0.00016831 |                   0.00078883 |                         8 |
|  8 | GLTSCR1L |                 0.00018744 |                   0.00087808 |                         9 |
|  9 | SLC2A1   |                 0.00018895 |                   0.00088384 |                        10 |


# ---
article_title: {{MND1}} and {{PSMC3IP}} Control {{PARP}} Inhibitor Sensitivity in Mitotic Cells

dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1J_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['Z-score/Z-ratios']

conditions: ['Olaparib']

sample:
|    | GeneID   |   CRISPRn Viability Z-score |   CRISPRn Olaparib Effect Z-score |   CRISPRn Olaparib sgRNA Z-count |
|---:|:---------|----------------------------:|----------------------------------:|---------------------------------:|
|  0 | A1BG     |                   -0.67095  |                         0.0944228 |                                0 |
|  1 | A1CF     |                    0.758922 |                        -4.03972   |                                1 |
|  2 | A2M      |                   -0.213711 |                         0.169044  |                                0 |
|  3 | A2ML1    |                    0.529894 |                         1.83532   |                                0 |
|  4 | A3GALT2  |                   -1.09668  |                         0.895692  |                                0 |
|  5 | A4GALT   |                   -0.464787 |                         0.496808  |                                0 |
|  6 | A4GNT    |                   -0.13382  |                         0.256893  |                                0 |
|  7 | AAAS     |                    0.125844 |                         0.291341  |                                0 |
|  8 | AACS     |                   -0.851788 |                         1.11655   |                                0 |
|  9 | AADAC    |                   -0.210665 |                        -1.85882   |                                1 |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Vincristine_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Vincristine']

sample:
|    | id      |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:--------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | TP53    |     6 |     0.97305 |       0.97497 |         1 |      19574 |               1 |   0.93454 |    2.44e-12 |      2.42e-07 |   0.00045 |          1 |               5 |   0.93454 |
|  1 | TRMT112 |     6 |     1       |       1       |         1 |      20458 |               0 |   0.9559  |    3.77e-11 |      2.42e-07 |   0.00045 |          2 |               5 |   0.9559  |
|  2 | RRP9    |     6 |     1       |       1       |         1 |      20466 |               0 |   0.67563 |    4.62e-11 |      2.42e-07 |   0.00045 |          3 |               6 |   0.67563 |
|  3 | CRCP    |     6 |     0.99904 |       0.99904 |         1 |      20295 |               0 |   0.91053 |    1.72e-10 |      2.42e-07 |   0.00045 |          4 |               5 |   0.91053 |
|  4 | NOP16   |     6 |     1       |       1       |         1 |      20463 |               0 |   0.99505 |    1.82e-10 |      2.42e-07 |   0.00045 |          5 |               6 |   0.99505 |
|  5 | LIPT1   |     6 |     1       |       1       |         1 |      20465 |               0 |   0.46702 |    5.51e-10 |      2.42e-07 |   0.00045 |          6 |               6 |   0.46702 |
|  6 | MLST8   |     6 |     1       |       1       |         1 |      20464 |               0 |   0.6126  |    6.64e-10 |      2.42e-07 |   0.00045 |          7 |               6 |   0.6126  |
|  7 | TP53BP1 |     6 |     0.99026 |       0.99021 |         1 |      19960 |               0 |   0.62782 |    2.04e-09 |      2.42e-07 |   0.00045 |          8 |               5 |   0.62782 |
|  8 | USP36   |     6 |     1       |       1       |         1 |      20462 |               0 |   0.3788  |    6.78e-09 |      2.42e-07 |   0.00045 |          9 |               6 |   0.3788  |
|  9 | EIF3I   |     6 |     0.99874 |       0.99872 |         1 |      20275 |               0 |   0.89706 |    6.8e-09  |      2.42e-07 |   0.00045 |         10 |               5 |   0.89706 |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Daunorubicin_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Daunorubicin']

sample:
|    | id     |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:-------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | HINFP  |     6 |     1       |       1       |         1 |      20466 |               0 |   1.8969  |  5.9965e-16 |    2.4189e-07 |  0.000248 |          1 |               6 |   1.8969  |
|  1 | TOP2B  |     6 |     0.95021 |       0.95021 |         1 |      18965 |               0 |   2.4823  |  4.504e-13  |    2.4189e-07 |  0.000248 |          2 |               5 |   2.4823  |
|  2 | CSTF3  |     6 |     1       |       1       |         1 |      20465 |               0 |   1.1544  |  3.9642e-12 |    2.4189e-07 |  0.000248 |          3 |               6 |   1.1544  |
|  3 | CRCP   |     6 |     1       |       1       |         1 |      20463 |               0 |   1.6789  |  2.6661e-11 |    2.4189e-07 |  0.000248 |          4 |               6 |   1.6789  |
|  4 | IARS2  |     6 |     1       |       1       |         1 |      20464 |               0 |   1.5182  |  4.2423e-11 |    2.4189e-07 |  0.000248 |          5 |               6 |   1.5182  |
|  5 | VPRBP  |     6 |     0.99982 |       0.99983 |         1 |      20364 |               0 |   1.1907  |  1.1185e-10 |    2.4189e-07 |  0.000248 |          6 |               6 |   1.1907  |
|  6 | NOP16  |     6 |     1       |       1       |         1 |      20462 |               0 |   1.3466  |  1.2664e-10 |    2.4189e-07 |  0.000248 |          7 |               6 |   1.3466  |
|  7 | CSNK2B |     6 |     0.90558 |       0.90562 |         1 |      18017 |               0 |   1.3907  |  1.3191e-10 |    2.4189e-07 |  0.000248 |          8 |               5 |   1.3907  |
|  8 | TP53   |     6 |     0.86402 |       0.86389 |         1 |      17166 |               0 |   0.93208 |  6.3852e-10 |    2.4189e-07 |  0.000248 |          9 |               5 |   0.93208 |
|  9 | DHDDS  |     6 |     1       |       1       |         1 |      20454 |               0 |   0.94588 |  6.6039e-10 |    2.4189e-07 |  0.000248 |         10 |               6 |   0.94588 |


# ---
article_title: Loss of the Abasic Site Sensor {{HMCES}} Is Synthetic Lethal with the Activity of the {{APOBEC3A}} Cytosine Deaminase in Cancer Cells

dataset_filename: biaynaLossAbasicSite2021_33788831_pbio.3001176.s014_lfc_all_time_points_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['fold-change/log fold-change']

conditions: ['A549 TP53-/-', 'A549 TP53 WT', 'LXF289 DOX-IC25', 'LXF289 DOX-IC50']

sample:
|    | Unnamed: 0   | A549 TP53-/-   | Unnamed: 2   | Unnamed: 3   | A549 TP53 WT   | Unnamed: 5   | Unnamed: 6   | LXF289 DOX-IC25   | Unnamed: 8   | Unnamed: 9   | LXF289 DOX-IC50   | Unnamed: 11   | Unnamed: 12   |
|---:|:-------------|:---------------|:-------------|:-------------|:---------------|:-------------|:-------------|:------------------|:-------------|:-------------|:------------------|:--------------|:--------------|
|  0 | Gene         | time 9         | time 12      | time 15      | time 9         | time 12      | time 15      | time 5            | time 10      | time 15      | time 5            | time 10       | time 15       |
|  1 | A1BG         | 0.24838        | -0.26098     | 0.099242     | 0.35132        | -0.025552    | -0.01677     | -0.07502          | -0.073809    | -0.22727     | 0.051521          | 0.041967      | -0.069046     |
|  2 | A1CF         | -0.006276      | 0.0038166    | -0.002651    | -0.25883       | 0.07213      | 0.060406     | -0.0064161        | -0.084975    | 0.057137     | -0.1384           | -0.14075      | 0.1142        |
|  3 | A2M          | 0.12738        | -0.1451      | -0.12855     | 0.19655        | 0.072992     | -0.012698    | 0.084314          | 0.28355      | -0.013266    | -0.15669          | 0.0010223     | -0.51692      |
|  4 | A2ML1        | -0.0030069     | -0.039999    | -0.064889    | 0.13849        | 0.11629      | 0.07048      | -0.1096           | -0.022392    | 0.2551       | -0.19769          | 0.051254      | 0.24108       |
|  5 | A3GALT2      | 0.053172       | 0.16654      | 0.089998     | 0.23265        | 0.03627      | 0.29674      | -0.089243         | -0.023538    | -0.17148     | -0.1837           | -0.069794     | -0.13972      |
|  6 | A4GALT       | 0.12626        | 0.22337      | 0.41143      | -0.036876      | -0.31352     | -0.50471     | -0.22779          | 0.10711      | 0.15743      | -0.028476         | -0.01239      | -0.2166       |
|  7 | A4GNT        | -0.076745      | 0.14238      | -0.24434     | -0.073255      | 0.14613      | -0.052216    | 0.0078288         | 0.2541       | 0.28358      | -0.057654         | 0.10609       | 0.15515       |
|  8 | AAAS         | 0.079229       | 0.029306     | -0.095838    | 0.253          | -0.087688    | -0.14106     | -0.011844         | -0.34891     | 0.011074     | -0.095835         | -0.25964      | -0.14949      |
|  9 | AACS         | 0.25882        | 0.017451     | 0.054099     | 0.084077       | 0.34883      | -0.63024     | 0.024298          | -0.090265    | 0.23925      | 0.0038947         | -0.030683     | 0.0038235     |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_6398_DMSO_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['DMSO']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   474 |
|---:|:----------------------------|------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |  1175 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |    82 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |   573 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |  1212 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |  1374 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |   376 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |   214 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |   467 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |     8 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |   503 |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_6377_Library_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   1021 |
|---:|:----------------------------|-------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |   1485 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |     26 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |   1665 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |   1068 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |    802 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |    976 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |    608 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |    374 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |    533 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |    611 |


# ---
article_title: Pooled CRISPR Screening in Pancreatic Cancer Cells Implicates Co-Repressor Complexes as a Cause of Multiple Drug Resistance via Regulation of Epithelial-to-Mesenchymal Transition

dataset_filename: ramakerPooledCRISPRScreening2021_34049503_12885_2021_8388_MOESM1_ESM_TableS3A_t1

header_row: 2

sgRNA_sequence: None

gene_identifier: RefSeq ID

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['fold-change/log fold-change']

conditions: ['Gemcitabine', 'Irinotecan', 'Oxaliplatin', '5FU']

sample:
|    | Supplemental Table S3A. Activation (SAM) Screen results. L2FC Sum for individual drug-cell line combinations and combined   | Unnamed: 1   | Unnamed: 2        | Unnamed: 3   | Unnamed: 4   | Unnamed: 5              | Unnamed: 6       | Unnamed: 7         | Unnamed: 8             | Unnamed: 9      | Unnamed: 10       | Unnamed: 11             | Unnamed: 12      | Unnamed: 13        | Unnamed: 14     | Unnamed: 15   | Unnamed: 16   | Unnamed: 17             | Unnamed: 18      | Unnamed: 19        | Unnamed: 20            | Unnamed: 21     | Unnamed: 22       | Unnamed: 23             | Unnamed: 24      | Unnamed: 25        | Unnamed: 26     | Unnamed: 27   | Unnamed: 28   |
|---:|:----------------------------------------------------------------------------------------------------------------------------|:-------------|:------------------|:-------------|:-------------|:------------------------|:-----------------|:-------------------|:-----------------------|:----------------|:------------------|:------------------------|:-----------------|:-------------------|:----------------|:--------------|:--------------|:------------------------|:-----------------|:-------------------|:-----------------------|:----------------|:------------------|:------------------------|:-----------------|:-------------------|:----------------|:--------------|:--------------|
|  0 | P1= Panc1,B3=BxPC3                                                                                                          | nan          | nan               | nan          | nan          | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan           | nan           | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan           | nan           |
|  1 | For each gene in each screen we calculated a L2FC sum, a p-value and FDR                                                    | nan          | nan               | nan          | nan          | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan           | nan           | nan                     | nan              | nan                | nan                    | nan             | nan               | nan                     | nan              | nan                | nan             | nan           | nan           |
|  2 | Gene_ID                                                                                                                     | nan          | Combined_L2FC_Sum | Combined_P   | Combined_FDR | P1_Gemcitabine_L2FC_Sum | P1_Gemcitabine_P | P1_Gemcitabine_FDR | P1_Irinotecan_L2FC_Sum | P1_Irinotecan_P | P1_Irinotecan_FDR | P1_Oxaliplatin_L2FC_Sum | P1_Oxaliplatin_P | P1_Oxaliplatin_FDR | P1_5FU_L2FC_Sum | P1_5FU_P      | P1_5FU_FDR    | B3_Gemcitabine_L2FC_Sum | B3_Gemcitabine_P | B3_Gemcitabine_FDR | B3_Irinotecan_L2FC_Sum | B3_Irinotecan_P | B3_Irinotecan_FDR | B3_Oxaliplatin_L2FC_Sum | B3_Oxaliplatin_P | B3_Oxaliplatin_FDR | B3_5FU_L2FC_Sum | B3_5FU_P      | B3_5FU_FDR    |
|  3 | NM_000014                                                                                                                   | A2M          | -0.18554717       | 0.130953394  | 0.350599116  | 1.30817191              | 0.02168621       | 0.129813642        | -0.093118004           | 0.050648391     | 0.16322417        | 1.565240464             | 0.004759342      | 0.048339743        | 0.31002847      | 0.171653606   | 0.394650874   | -0.087805172            | 0.594800367      | 0.890740617        | -0.170216254           | 0.62710939      | 0.916172815       | -0.147466101            | 0.669739761      | 0.911610932        | -0.336701153    | 0.744820208   | 0.948290624   |
|  4 | NM_000015                                                                                                                   | NAT2         | -0.615084931      | 0.226631071  | 0.380142049  | 0.404445504             | 0.218215912      | 0.461666725        | -2.04439538            | 0.231691291     | 0.428193874       | -1.161697916            | 0.827891999      | 0.944957403        | 0.341308066     | 0.160581291   | 0.379965492   | 0.163152039             | 0.252576833      | 0.693515925        | 0.306888492            | 0.131873646     | 0.603273276       | 0.520695716             | 0.021178454      | 0.285921004        | -0.246447184    | 0.665934009   | 0.921367956   |
|  5 | NM_000017                                                                                                                   | ACADS        | -0.098666736      | 0.114730206  | 0.350599116  | 1.349700228             | 0.018844047      | 0.119126314        | -0.592403999           | 0.081263289     | 0.220288023       | 1.179154702             | 0.016340687      | 0.102186584        | -0.345339518    | 0.479544188   | 0.693959826   | -0.165334204            | 0.692739101      | 0.92914601         | -0.499325813           | 0.883811164     | 0.98533932        | 0.109755175             | 0.320029306      | 0.72789716         | 0.1602379       | 0.245797754   | 0.684280521   |
|  6 | NM_000018                                                                                                                   | ACADVL       | -0.32636044       | 0.159638139  | 0.350599116  | -0.262644238            | 0.517900163      | 0.72989745         | -2.127690308           | 0.24271302      | 0.441365398       | 0.439541554             | 0.122595823      | 0.342688692        | 0.645351231     | 0.07637892    | 0.247028188   | -0.096276209            | 0.606091263      | 0.895640935        | 0.120764062            | 0.295917814     | 0.762533957       | 0.128628439             | 0.295466635      | 0.70867606         | 0.112005262     | 0.291027849   | 0.722306989   |
|  7 | NM_000019                                                                                                                   | ACAT1        | -0.725333818      | 0.254138716  | 0.399478635  | -0.169671931            | 0.47495659       | 0.697794863        | -2.202770499           | 0.252692122     | 0.452401402       | -0.838159859            | 0.710191083      | 0.885133611        | 0.309267018     | 0.171942817   | 0.394882426   | 0.042204207             | 0.412285297      | 0.801898115        | 0.377173277            | 0.090371325     | 0.537218855       | 0.584798233             | 0.01146539       | 0.234084054        | -0.032607272    | 0.441387353   | 0.820930116   |
|  8 | NM_000020                                                                                                                   | ACVRL1       | 0.024216187       | 0.093877587  | 0.350599116  | -2.427987483            | 0.960279087      | 0.984455555        | 6.74630897             | 2.28e-06        | 0.000329514       | -0.256757607            | 0.417704253      | 0.681421572        | 0.193727563     | 0.217458582   | 0.453125559   | -0.033466217            | 0.519405286      | 0.855408717        | -0.14850765            | 0.603527978     | 0.909010688       | -0.126252199            | 0.64368972       | 0.901530275        | 0.405090813     | 0.08247522    | 0.468419808   |
|  9 | NM_000021                                                                                                                   | PSEN1        | -0.110484192      | 0.116867676  | 0.350599116  | 0.328251514             | 0.248306685      | 0.492817485        | -0.86647071            | 0.102642842     | 0.256236051       | 0.805886007             | 0.049180567      | 0.201025893        | 0.464599922     | 0.121267711   | 0.325521892   | -0.094310772            | 0.603531993      | 0.894589312        | 0.38642245             | 0.085724494     | 0.525854217       | -0.17994996             | 0.707301358      | 0.924954734        | -0.554098487    | 0.881931648   | 0.985445813   |


# ---
article_title: {{CCNE1}} Amplification Is Synthetic Lethal with {{PKMYT1}} Kinase Inhibition

dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts RP-6306 resistance_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['RP-6306 resistance']

sample:
|    | Read counts for FT282-hTERT p53-R175H CCNE1-high RP-6306 Resistance Screen for Clones 3 and 4   | Unnamed: 1   | Unnamed: 2                    | Unnamed: 3                    | Unnamed: 4                    | Unnamed: 5                    | Unnamed: 6                         | Unnamed: 7                         | Unnamed: 8                         | Unnamed: 9                         |
|---:|:------------------------------------------------------------------------------------------------|:-------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:-----------------------------------|:-----------------------------------|:-----------------------------------|:-----------------------------------|
|  0 | sgRNA                                                                                           | Gene         | TKOv3_FT282.CCNE1oe_C3_NT_T6A | TKOv3_FT282.CCNE1oe_C3_NT_T6B | TKOv3_FT282.CCNE1oe_C4_NT_T6A | TKOv3_FT282.CCNE1oe_C4_NT_T6B | TKOv3_FT282.CCNE1oe_C3_RP6306_T21A | TKOv3_FT282.CCNE1oe_C3_RP6306_T21B | TKOv3_FT282.CCNE1oe_C4_RP6306_T21A | TKOv3_FT282.CCNE1oe_C4_RP6306_T21B |
|  1 | chr11:134201957-134201976_GLB1L2_+                                                              | GLB1L2       | 1737                          | 1521                          | 1549                          | 1406                          | 708                                | 1100                               | 1178                               | 843                                |
|  2 | chr6:26022026-26022045_HIST1H4A_-                                                               | HIST1H4A     | 2052                          | 1379                          | 2108                          | 1631                          | 1735                               | 975                                | 1723                               | 853                                |
|  3 | chr12:120436384-120436403_CCDC64_+                                                              | CCDC64       | 425                           | 391                           | 496                           | 476                           | 146                                | 506                                | 337                                | 92                                 |
|  4 | chr8:73979678-73979697_SBSPON_-                                                                 | SBSPON       | 894                           | 800                           | 1048                          | 780                           | 956                                | 1193                               | 1100                               | 1100                               |
|  5 | chr13:28636131-28636150_FLT3_-                                                                  | FLT3         | 1329                          | 1147                          | 1677                          | 1370                          | 795                                | 949                                | 1549                               | 1848                               |
|  6 | chr16:27549581-27549600_GTF3C1_+                                                                | GTF3C1       | 669                           | 546                           | 753                           | 660                           | 753                                | 438                                | 623                                | 342                                |
|  7 | chr7:87082265-87082284_ABCB4_+                                                                  | ABCB4        | 2140                          | 2598                          | 2661                          | 2156                          | 1554                               | 1012                               | 2428                               | 3655                               |
|  8 | chr5:37834809-37834828_GDNF_+                                                                   | GDNF         | 962                           | 764                           | 796                           | 1015                          | 598                                | 2018                               | 1414                               | 2447                               |
|  9 | chr9:97535429-97535448_C9orf3_+                                                                 | C9orf3       | 960                           | 719                           | 877                           | 960                           | 323                                | 551                                | 1479                               | 550                                |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_6397_Puromycin_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   373 |
|---:|:----------------------------|------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |   814 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |    10 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |   344 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |   880 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |  1049 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |    84 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |   220 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |   232 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |     5 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |   283 |


# ---
article_title: Identifying CDC7 as a Synergistic Target of Chemotherapy in Resistant Small-Cell Lung Cancer via CRISPR/Cas9 Screening

dataset_filename: dengIdentifyingCDC7Synergistic2023_36725843_41420_2023_1315_MOESM8_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['d7_chemo', 'd7_control']

sample:
|    | id       |   d7_chemo__d7_control__num |   d7_chemo__d7_control__neg_p-value |   d7_chemo__d7_control__neg_fdr |   d7_chemo__d7_control__neg_rank |   d7_chemo__d7_control__neg_goodsgrna |   d7_chemo__d7_control__pos_p-value |   d7_chemo__d7_control__pos_fdr |   d7_chemo__d7_control__pos_rank |   d7_chemo__d7_control__pos_goodsgrna |   d7_chemo__d7_control__log2FoldChange |
|---:|:---------|----------------------------:|------------------------------------:|--------------------------------:|---------------------------------:|--------------------------------------:|------------------------------------:|--------------------------------:|---------------------------------:|--------------------------------------:|---------------------------------------:|
|  0 | TRMT12   |                           3 |                             0.29732 |                        0.996322 |                             6050 |                                     1 |                          5.3276e-06 |                        0.113861 |                                1 |                                     2 |                                 9.1497 |
|  1 | C11orf16 |                           1 |                             0.99868 |                        0.999829 |                            21344 |                                     0 |                          0.0013627  |                        0.956958 |                               68 |                                     1 |                                 9.0236 |
|  2 | FAM208B  |                           1 |                             0.99693 |                        0.999632 |                            21312 |                                     0 |                          0.0030884  |                        0.956958 |                              152 |                                     1 |                                 8.8196 |
|  3 | APOBR    |                           1 |                             0.996   |                        0.999471 |                            21297 |                                     0 |                          0.0039862  |                        0.956958 |                              196 |                                     1 |                                 8.7494 |
|  4 | WFDC13   |                           1 |                             0.99597 |                        0.999471 |                            21296 |                                     0 |                          0.004027   |                        0.956958 |                              197 |                                     1 |                                 8.7472 |
|  5 | ARHGDIG  |                           1 |                             0.99563 |                        0.999463 |                            21285 |                                     0 |                          0.0043387  |                        0.956958 |                              219 |                                     1 |                                 8.7226 |
|  6 | MYH9     |                           2 |                             0.99986 |                        0.999904 |                            21369 |                                     0 |                          0.00028839 |                        0.880481 |                               10 |                                     2 |                                 8.7127 |
|  7 | MISP     |                           1 |                             0.99513 |                        0.999463 |                            21276 |                                     0 |                          0.0048349  |                        0.956958 |                              245 |                                     1 |                                 8.664  |
|  8 | RBMY1J   |                           1 |                             0.99492 |                        0.999463 |                            21271 |                                     0 |                          0.0050448  |                        0.966452 |                              256 |                                     1 |                                 8.6415 |
|  9 | AMH      |                           3 |                             0.36741 |                        0.996322 |                             7300 |                                     1 |                          0.00026615 |                        0.880481 |                                4 |                                     2 |                                 8.5893 |


# ---
article_title: {{MND1}} and {{PSMC3IP}} Control {{PARP}} Inhibitor Sensitivity in Mitotic Cells

dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1C_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Z-score/Z-ratios']

conditions: ['olap screen', 'talaz screen']

sample:
|    | GeneID   | Ensembl gene ID   | sgRNA                 |   CRISPRi Viability Z score (olap screen) |   CRISPRi Drug Effect Z score (olap screen) |   CRISPRi Viability Z score (talaz screen) |   CRISPRi Drug Effect Z score (talaz screen) | Essentiality   |
|---:|:---------|:------------------|:----------------------|------------------------------------------:|--------------------------------------------:|-------------------------------------------:|---------------------------------------------:|:---------------|
|  0 | A1BG     | ENSG00000121410   | A1BG_+_58864367.23-P2 |                                  0.12549  |                                   0.470385  |                                   0.317818 |                                   -0.551194  | non-essential  |
|  1 | A1BG     | ENSG00000121410   | A1BG_-_58864840.23-P2 |                                  0.956319 |                                  -0.0779968 |                                   1.3407   |                                   -0.413111  | non-essential  |
|  2 | A1BG     | ENSG00000121410   | A1BG_-_58864822.23-P2 |                                 -0.451518 |                                   0.449761  |                                  -0.524998 |                                    0.637648  | non-essential  |
|  3 | A1BG     | ENSG00000121410   | A1BG_+_58858549.23-P1 |                                  0.805702 |                                  -0.589471  |                                   0.717305 |                                    0.0119589 | non-essential  |
|  4 | A1BG     | ENSG00000121410   | A1BG_-_58864360.23-P2 |                                 -1.05614  |                                   0.0701778 |                                  -1.2004   |                                   -0.628214  | non-essential  |
|  5 | A1BG     | ENSG00000121410   | A1BG_-_58858617.23-P1 |                                  1.21009  |                                  -1.32642   |                                   1.7404   |                                   -0.649347  | non-essential  |
|  6 | A1BG     | ENSG00000121410   | A1BG_-_58858788.23-P1 |                                 -0.3081   |                                   0.751719  |                                  -0.567585 |                                    2.04275   | non-essential  |
|  7 | A1BG     | ENSG00000121410   | A1BG_-_58858630.23-P1 |                                 -0.69372  |                                  -0.844938  |                                  -1.33596  |                                    0.804068  | non-essential  |
|  8 | A1BG     | ENSG00000121410   | A1BG_+_58858964.23-P1 |                                  0.515829 |                                   0.894974  |                                   0.499489 |                                    0.317354  | non-essential  |
|  9 | A1BG     | ENSG00000121410   | A1BG_-_58864354.23-P2 |                                  1.41769  |                                  -1.90994   |                                   0.991076 |                                   -0.620136  | non-essential  |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_L-asparaginase_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['L-asparaginase']

sample:
|    | id     |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:-------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | DHODH  |     6 |     0.98227 |       0.98211 |         1 |      19328 |               0 |    2.9941 |    1.09e-11 |      2.42e-07 |   0.00016 |          1 |               5 |    2.9941 |
|  1 | G6PD   |     6 |     1       |       1       |         1 |      20466 |               0 |    1.7795 |    1.37e-11 |      2.42e-07 |   0.00016 |          2 |               6 |    1.7795 |
|  2 | CAD    |     6 |     1       |       1       |         1 |      20465 |               0 |    1.8716 |    1.42e-11 |      2.42e-07 |   0.00016 |          3 |               6 |    1.8716 |
|  3 | PMPCA  |     6 |     1       |       1       |         1 |      20464 |               0 |    1.4795 |    2.08e-11 |      2.42e-07 |   0.00016 |          4 |               6 |    1.4795 |
|  4 | IMPDH2 |     6 |     1       |       1       |         1 |      20463 |               0 |    1.6361 |    6.2e-11  |      2.42e-07 |   0.00016 |          5 |               6 |    1.6361 |
|  5 | NAA10  |     6 |     0.99587 |       0.99585 |         1 |      19801 |               0 |    2.9033 |    6.45e-11 |      2.42e-07 |   0.00016 |          6 |               5 |    2.9033 |
|  6 | GMPPB  |     6 |     1       |       1       |         1 |      20462 |               0 |    2.1461 |    6.46e-11 |      2.42e-07 |   0.00016 |          7 |               6 |    2.1461 |
|  7 | COQ2   |     6 |     1       |       1       |         1 |      20461 |               0 |    1.642  |    2.75e-10 |      2.42e-07 |   0.00016 |          8 |               6 |    1.642  |
|  8 | EIF2S3 |     6 |     0.99999 |       0.99999 |         1 |      20315 |               0 |    2.5704 |    2.95e-10 |      2.42e-07 |   0.00016 |          9 |               6 |    2.5704 |
|  9 | FPGS   |     6 |     1       |       1       |         1 |      20460 |               0 |    1.7745 |    3.82e-10 |      2.42e-07 |   0.00016 |         10 |               6 |    1.7745 |


# ---
article_title: A Genome-Scale CRISPR/Cas9 Knockout Screening Reveals SH3D21 as a Sensitizer for Gemcitabine

dataset_filename: masoudiGenomescaleCRISPRCas92019_31844142_41598_2019_55893_MOESM2_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Gemcitabine']

sample:
|    | Gene                                  |   # Hairpins |      NES |   Gene rank |   p-value |   p-value rank |
|---:|:--------------------------------------|-------------:|---------:|------------:|----------:|---------------:|
|  0 | SH3D21                                |            6 | 0.003022 |           1 |    0.0001 |              2 |
|  1 | DNAJB12                               |            6 | 0.003166 |           2 |    0.0001 |              3 |
|  2 | hsa-mir-4443                          |            4 | 0.003691 |           3 |    0.0001 |              1 |
|  3 | COA3                                  |            6 | 0.005052 |           4 |    0.0003 |              5 |
|  4 | ATAD2                                 |            6 | 0.007434 |           5 |    0.0004 |              6 |
|  5 | UBA1                                  |            6 | 0.007651 |           6 |    0.0005 |              8 |
|  6 | NonTargetingControlGuideForHuman_0542 |            2 | 0.008692 |           7 |    0.0002 |              4 |
|  7 | COPA                                  |            6 | 0.008835 |           8 |    0.0005 |              9 |
|  8 | NOL11                                 |            6 | 0.009272 |           9 |    0.0007 |             13 |
|  9 | RAB17                                 |            6 | 0.009596 |          10 |    0.0007 |             12 |


# ---
article_title: Ribosomal Protein S11 Influences Glioma Response to TOP2 Poisons

dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Eto R1', 'Eto R2', 'DMSO']

sample:
|    | genes   |   Eto R1 |   Eto R2 |   DMSO |
|---:|:--------|---------:|---------:|-------:|
|  0 | A1BG    |     2032 |      286 |    677 |
|  1 | A1BG    |      934 |      497 |   1248 |
|  2 | A1BG    |        0 |        8 |     88 |
|  3 | A1BG    |     2544 |      726 |    842 |
|  4 | A1CF    |     1191 |      236 |    972 |
|  5 | A1CF    |      882 |      343 |    749 |
|  6 | A1CF    |     1188 |      399 |    509 |
|  7 | A1CF    |     1118 |      209 |    145 |
|  8 | A2M     |      954 |      231 |    865 |
|  9 | A2M     |     2044 |      482 |    677 |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_6396_Puromycin_Day0_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Puromycin']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   608 |
|---:|:----------------------------|------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |  1154 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |    31 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |  1008 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |   852 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |  1040 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |   352 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |   307 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |   406 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |    97 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |   283 |


# ---
article_title: Ribosomal Protein S11 Influences Glioma Response to TOP2 Poisons

dataset_filename: awahRibosomalProteinS112020_32528131_NIHMS1599130-supplement-1599130_SuppList_Data_6_CRISPR Read count Awah CU et.al_t3

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Library R1', 'Library R2']

sample:
|    | sgRNA                      |   Library R1 |   Library R2 |
|---:|:---------------------------|-------------:|-------------:|
|  0 | A1BG_CATCTTCTTTCACCTGAACG  |         1373 |         1004 |
|  1 | A1BG_CTCCGGGGAGAACTCCGGCG  |         1188 |         1496 |
|  2 | A1BG_TCTCCATGGTGCATCAGCAC  |           21 |           14 |
|  3 | A1BG_TGGAAGTCCACTCCACTCAG  |         2090 |         1652 |
|  4 | A1CF_ACAGGAAGAATTCAGTTATG  |         1101 |         1053 |
|  5 | A1CF_AGTTATGTTAGGTATACCCG  |          928 |          749 |
|  6 | A1CF_CTTCATTTCCCAGCCACCAA  |         1244 |          953 |
|  7 | A1CF_GAATTCAACAATATCAAACC  |          948 |          586 |
|  8 | A2ML1_AAATACTTACTGGTATCGAG |          878 |          405 |
|  9 | A2ML1_ACCCTGGTTACTGATAACAA |          374 |          483 |


# ---
article_title: {{CCNE1}} Amplification Is Synthetic Lethal with {{PKMYT1}} Kinase Inhibition

dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts TKOv2 CCNE1 SL_t1

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['TKOv2', 'RPE1.hTERT.Cas9.dTP53_CCNE1oe', 'Clone#2_T0', 'Clone#2_T15_repA', 'Clone#2_T15_repB', 'Clone#2_T18_repA', 'Clone#2_T18_repB', 'Clone#21_T0', 'Clone#21_T15_repA', 'Clone#21_T15_repB', 'Clone#21_T18_repA', 'Clone#21_T18_repB', 'RPE1.hTERT.Cas9.dTP53_T0', 'RPE1.hTERT.Cas9.dTP53_WT_T18_repA', 'RPE1.hTERT.Cas9.dTP53_T18_repB']

sample:
|    | Read Counts for screens using TKOv2   | Unnamed: 1   | Unnamed: 2                                     | Unnamed: 3                                           | Unnamed: 4                                           | Unnamed: 5                                           | Unnamed: 6                                           | Unnamed: 7                                      | Unnamed: 8                                            | Unnamed: 9                                            | Unnamed: 10                                           | Unnamed: 11                                           | Unnamed: 12                    | Unnamed: 13                             | Unnamed: 14                          |
|---:|:--------------------------------------|:-------------|:-----------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:-----------------------------------------------------|:------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:------------------------------------------------------|:-------------------------------|:----------------------------------------|:-------------------------------------|
|  0 | sgRNA                                 | Gene         | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T15_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T15_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T15_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T15_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#21_T18_repB | TKOv2_RPE1.hTERT.Cas9.dTP53_T0 | TKOv2_RPE1.hTERT.Cas9.dTP53_WT_T18_repA | TKOv2_RPE1.hTERT.Cas9.dTP53_T18_repB |
|  1 | chr1:100111914-100111933_PALMD_-      | PALMD        | 516                                            | 242                                                  | 166                                                  | 325                                                  | 283                                                  | 1106                                            | 220                                                   | 521                                                   | 305                                                   | 117                                                   | 1340                           | 690                                     | 676                                  |
|  2 | chr1:100133222-100133241_PALMD_+      | PALMD        | 101                                            | 35                                                   | 31                                                   | 52                                                   | 37                                                   | 129                                             | 72                                                    | 99                                                    | 50                                                    | 35                                                    | 167                            | 117                                     | 62                                   |
|  3 | chr1:100152265-100152284_PALMD_+      | PALMD        | 808                                            | 290                                                  | 134                                                  | 370                                                  | 255                                                  | 703                                             | 194                                                   | 902                                                   | 218                                                   | 179                                                   | 1316                           | 705                                     | 586                                  |
|  4 | chr1:100154717-100154736_PALMD_+      | PALMD        | 425                                            | 239                                                  | 133                                                  | 170                                                  | 201                                                  | 476                                             | 122                                                   | 453                                                   | 145                                                   | 153                                                   | 439                            | 421                                     | 277                                  |
|  5 | chr1:100203758-100203777_FRRS1_-      | FRRS1        | 811                                            | 148                                                  | 135                                                  | 256                                                  | 197                                                  | 514                                             | 171                                                   | 513                                                   | 129                                                   | 108                                                   | 876                            | 452                                     | 395                                  |
|  6 | chr1:100206368-100206387_FRRS1_+      | FRRS1        | 915                                            | 136                                                  | 309                                                  | 247                                                  | 353                                                  | 1256                                            | 194                                                   | 821                                                   | 168                                                   | 220                                                   | 1480                           | 1287                                    | 628                                  |
|  7 | chr1:100212897-100212916_FRRS1_+      | FRRS1        | 176                                            | 57                                                   | 58                                                   | 105                                                  | 40                                                   | 249                                             | 28                                                    | 197                                                   | 50                                                    | 77                                                    | 580                            | 147                                     | 253                                  |
|  8 | chr1:100214146-100214165_FRRS1_-      | FRRS1        | 604                                            | 146                                                  | 199                                                  | 114                                                  | 199                                                  | 509                                             | 80                                                    | 369                                                   | 100                                                   | 95                                                    | 1226                           | 422                                     | 235                                  |
|  9 | chr1:100316659-100316678_AGL_+        | AGL          | 1389                                           | 351                                                  | 513                                                  | 532                                                  | 449                                                  | 1989                                            | 367                                                   | 1216                                                  | 319                                                   | 262                                                   | 2876                           | 1151                                    | 735                                  |


# ---
article_title: Genome-Scale CRISPR Screen Reveals Neddylation to Contribute to Cisplatin Resistance of Testicular Germ Cell Tumours

dataset_filename: funkeGenomescaleCRISPRScreen2023_37024667_41416_2023_2247_MOESM3_ESM_7_NGS_Cand_list_t2

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Cisplatin']

sample:
|    |   Unnamed: 0 | Gene Name   | Refseq ID    | Target sequence sgRNA   |   Read counts |
|---:|-------------:|:------------|:-------------|:------------------------|--------------:|
|  0 |          nan | CYP2R1      | NM_024514    | AGAAGGTAGTTGTCCCGAAG    |       4248699 |
|  1 |          nan | LEKR1       | NM_001193283 | GGAGCGTTGACGACAGACGT    |       4084952 |
|  2 |          nan | CCL1        | NM_002981    | GAGTTTACCACATGGTCGCT    |       2230750 |
|  3 |          nan | LHB         | NM_000894    | CATGTGCACCTCTCGCCCCC    |       1223512 |
|  4 |          nan | WBP2NL      | NM_152613    | AGCGCGAAGGGGTGAGTGAT    |       1066313 |
|  5 |          nan | OR2W1       | NM_030903    | AAGGAAGAAGAAGGAAGAAG    |        331952 |
|  6 |          nan | PLEKHF1     | NM_024310    | CCCGAGCTTGTGCGCGCTGC    |        307083 |
|  7 |          nan | NAE1        | NM_003905    | GTCTGCCTTTCTGGCCGATG    |        179370 |
|  8 |          nan | CLDND2      | NM_152353    | TTAGATTCAGATCTCAGCTT    |        171947 |
|  9 |          nan | ANGPTL1     | NM_004673    | AGAGATTGTGTGCAGGTTGA    |        167831 |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Methotrexate_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Methotrexate']

sample:
|    | id            |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:--------------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | SLC19A1       |     6 |     0.97481 |       0.97475 |         1 |      19354 |               0 |    4.8868 |  6.637e-21  |    2.4189e-07 |  8e-05    |          1 |               5 |    4.8868 |
|  1 | NAA60         |     6 |     0.99972 |       0.99973 |         1 |      20257 |               0 |    4.4439 |  1.0139e-17 |    2.4189e-07 |  8e-05    |          2 |               6 |    4.4439 |
|  2 | PPP2CA        |     5 |     0.99927 |       0.99929 |         1 |      20200 |               0 |    3.2935 |  3.6789e-09 |    2.4189e-07 |  8e-05    |         45 |               5 |    3.2935 |
|  3 | CNOT4         |     6 |     1       |       1       |         1 |      20405 |               0 |    3.0542 |  3.6357e-16 |    2.4189e-07 |  8e-05    |          3 |               6 |    3.0542 |
|  4 | HNRNPM        |     6 |     0.98809 |       0.98801 |         1 |      19724 |               0 |    2.9752 |  8.8629e-13 |    2.4189e-07 |  8e-05    |          8 |               5 |    2.9752 |
|  5 | CAD           |     6 |     1       |       1       |         1 |      20463 |               0 |    2.9751 |  1.8228e-12 |    2.4189e-07 |  8e-05    |          9 |               6 |    2.9751 |
|  6 | hsa-mir-548ak |     1 |     0.99557 |       0.99561 |         1 |      19992 |               0 |    2.9244 |  0.0044338  |    0.017554   |  0.337848 |       1063 |               1 |    2.9244 |
|  7 | DHODH         |     6 |     0.99471 |       0.99473 |         1 |      19956 |               0 |    2.9139 |  6.6238e-11 |    2.4189e-07 |  8e-05    |         22 |               5 |    2.9139 |
|  8 | hsa-mir-3669  |     2 |     0.99864 |       0.99868 |         1 |      20156 |               0 |    2.8685 |  0.0010728  |    0.0047761  |  0.134084 |        729 |               2 |    2.8685 |
|  9 | NAA30         |     6 |     0.99991 |       0.99991 |         1 |      20319 |               0 |    2.8596 |  4.6508e-14 |    2.4189e-07 |  8e-05    |          5 |               6 |    2.8596 |


# ---
article_title: Pooled CRISPR Screening in Pancreatic Cancer Cells Implicates Co-Repressor Complexes as a Cause of Multiple Drug Resistance via Regulation of Epithelial-to-Mesenchymal Transition

dataset_filename: ramakerPooledCRISPRScreening2021_34049503_12885_2021_8388_MOESM1_ESM_TableS3B_t1

header_row: 2

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['fold-change/log fold-change']

conditions: ['Gemcitabine', 'Oxaliplatin', 'Irinotecan', '5FU']

sample:
|    | Supplemental Table S3B. Knockout (Gecko) Screen results. L2FC Sum for individual drug-cell line combinations    | Unnamed: 1         | Unnamed: 2         | Unnamed: 3        | Unnamed: 4              | Unnamed: 5         | Unnamed: 6         | Unnamed: 7              | Unnamed: 8         | Unnamed: 9         | Unnamed: 10            | Unnamed: 11        | Unnamed: 12       | Unnamed: 13        | Unnamed: 14          | Unnamed: 15       | Unnamed: 16             | Unnamed: 17       | Unnamed: 18        | Unnamed: 19             | Unnamed: 20       | Unnamed: 21        | Unnamed: 22            | Unnamed: 23        | Unnamed: 24       | Unnamed: 25        | Unnamed: 26        | Unnamed: 27       |
|---:|:----------------------------------------------------------------------------------------------------------------|:-------------------|:-------------------|:------------------|:------------------------|:-------------------|:-------------------|:------------------------|:-------------------|:-------------------|:-----------------------|:-------------------|:------------------|:-------------------|:---------------------|:------------------|:------------------------|:------------------|:-------------------|:------------------------|:------------------|:-------------------|:-----------------------|:-------------------|:------------------|:-------------------|:-------------------|:------------------|
|  0 | P1 = Panc1: B3=BxPC3                                                                                            | nan                | nan                | nan               | nan                     | nan                | nan                | nan                     | nan                | nan                | nan                    | nan                | nan               | nan                | nan                  | nan               | nan                     | nan               | nan                | nan                     | nan               | nan                | nan                    | nan                | nan               | nan                | nan                | nan               |
|  1 | For each gene in each screen we calculated a L2FC sum, a p-value and FDR                                        | nan                | nan                | nan               | nan                     | nan                | nan                | nan                     | nan                | nan                | nan                    | nan                | nan               | nan                | nan                  | nan               | nan                     | nan               | nan                | nan                     | nan               | nan                | nan                    | nan                | nan               | nan                | nan                | nan               |
|  2 | nan                                                                                                             | Combined_L2FC_Sum  | Combined_P         | Combined_FDR      | P1_Gemcitabine_L2FC_Sum | P1_Gemcitabine_P   | P1_Gemcitabine_FDR | P1_Oxaliplatin_L2FC_Sum | P1_Oxaliplatin_P   | P1_Oxaliplatin_FDR | P1_Irinotecan_L2FC_Sum | P1_Irinotecan_P    | P1_Irinotecan_FDR | P1_5FU_L2FC_Sum    | P1_5FU_P             | P1_5FU_FDR        | B3_Gemcitabine_L2FC_Sum | B3_Gemcitabine_P  | B3_Gemcitabine_FDR | B3_Oxaliplatin_L2FC_Sum | B3_Oxaliplatin_P  | B3_Oxaliplatin_FDR | B3_Irinotecan_L2FC_Sum | B3_Irinotecan_P    | B3_Irinotecan_FDR | B3_5FU_L2FC_Sum    | B3_5FU_P           | B3_5FU_FDR        |
|  3 | A1BG                                                                                                            | -0.904586216999054 | 0.314775208036283  | 0.701825580880729 | -0.904586216999054      | 0.705751098984387  | 0.999965452086164  | -0.392715133913581      | 0.72360982484397   | 0.970318452715167  | -1.4862571781002       | 0.28926700246282   | 0.999898342971155 | -2.04481219892761  | 0.618499289440127    | 0.99991789536423  | -0.378696386660878      | 0.747680849435469 | 0.993359561111788  | -0.00544005991251131    | 0.394241943897826 | 0.962196700272364  | -0.607283237222274     | 0.781932864584173  | 0.999628192829801 | -1.11520453166738  | 0.950561172993587  | 0.999687840103799 |
|  4 | A1CF                                                                                                            | 0.351737973373145  | 0.0193697286062979 | 0.650760987670613 | 0.351737973373145       | 0.0927310393613259 | 0.922223015712754  | 0.517615043804666       | 0.0444151399234951 | 0.618721214655904  | -0.107825095499609     | 0.0502095521796057 | 0.807361064319172 | 0.37677409411936   | 0.0518009297082954   | 0.689169112092335 | -0.684541032563738      | 0.903139463221877 | 0.997791930332779  | -0.362510690205391      | 0.673465770503558 | 0.981022096637038  | -0.378322136192466     | 0.641481478495057  | 0.999628192829801 | 0.510627260795835  | 0.0663233750561576 | 0.885783329272473 |
|  5 | A2M                                                                                                             | -2.30096680456106  | 0.754005655042413  | 0.88722863851172  | -2.30096680456106       | 0.947616593401041  | 0.999965452086164  | -0.419329703150225      | 0.743358163881619  | 0.974014852908947  | -1.58911651319461      | 0.317460927776845  | 0.999898342971155 | -0.600966589391017 | 0.28382254856844     | 0.927929483391031 | 0.190931131351119       | 0.231185979158862 | 0.95314044554226   | 0.295867372305665       | 0.180498180214185 | 0.957064868238283  | 0.396586814896698      | 0.0951293360641197 | 0.991619292435505 | 0.336148755140554  | 0.133589593267939  | 0.925925912587347 |
|  6 | A2ML1                                                                                                           | -1.9759314381047   | 0.676834867928176  | 0.851804423894426 | 0.442851932499182       | 0.0694194330756404 | 0.909520276361664  | -0.485279453354449      | 0.787588081336823  | 0.981000334589683  | -1.94732080198089      | 0.421180656036944  | 0.999898342971155 | -0.700232946824223 | 0.311817188453453    | 0.949261561256625 | -0.243154531677502      | 0.638762340942322 | 0.988923920805435  | -0.503046371740626      | 0.760450928221686 | 0.987683709242053  | -0.0404659417536793    | 0.366904400974052  | 0.991619292435505 | -0.129795195538754 | 0.472986939204748  | 0.989267682706392 |
|  7 | A3GALT2                                                                                                         | -0.820002762004142 | 0.283324069374965  | 0.690813793378523 | -0.820002762004142      | 0.671848820170785  | 0.999965452086164  | 0.111026558795745       | 0.260100664384941  | 0.851555112410188  | -4.50233933019628      | 0.945767369760943  | 0.999898342971155 | 1.84863278230968   | 0.000952819229919335 | 0.689169112092335 | -0.284740414302769      | 0.674853681034726 | 0.991320679217531  | -1.1465405039203        | 0.953036031154894 | 0.998162255674476  | -0.215192339078237     | 0.512693721879082  | 0.993282096401363 | -0.102168727299402 | 0.449133813216324  | 0.986372560408026 |
|  8 | A4GALT                                                                                                          | -0.44094605396728  | 0.155636514496081  | 0.658472760662647 | -0.44094605396728       | 0.483272194431812  | 0.999965452086164  | -0.239340864911618      | 0.593370243607812  | 0.944836859947227  | -1.95124727974868      | 0.422364509490053  | 0.999898342971155 | 0.033984583384969  | 0.110411475998078    | 0.736070079130959 | 0.0108939512195714      | 0.392217401987441 | 0.968386871502198  | -0.512376918742163      | 0.765486306791517 | 0.987683709242053  | -0.424077135816083     | 0.673517150736183  | 0.999628192829801 | 0.129286847058236  | 0.261175612596223  | 0.967679494098821 |
|  9 | A4GNT                                                                                                           | -0.941392782767956 | 0.328699340245052  | 0.707086673221779 | -0.941392782767956      | 0.719793845687434  | 0.999965452086164  | 0.00340518163832562     | 0.356760620092611  | 0.89066405225398   | -1.23742319608632      | 0.227079281755098  | 0.999898342971155 | -2.70896664336663  | 0.721284157048936    | 0.99991789536423  | -0.434111879562378      | 0.785470013171728 | 0.994972178394953  | 0.0130169457431446      | 0.3794580724168   | 0.962196700272364  | -1.02181254573655      | 0.918834362753794  | 0.999628192829801 | -1.25121690918341  | 0.965761496944825  | 0.999687840103799 |


# ---
article_title: Genome Scale CRISPR Cas9a Knockout Screen Reveals Genes That Control Glioblastoma Susceptibility to the Alkylating Agent Temozolomide

dataset_filename: awahGenomeScaleCRISPR2022_35990011_6399_Temozolomide_Day14_sheet_t1

header_row: None

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Temozolomide']

sample:
|    | A1BG_CATCTTCTTTCACCTGAACG   |   311 |
|---:|:----------------------------|------:|
|  0 | A1BG_CTCCGGGGAGAACTCCGGCG   |   355 |
|  1 | A1BG_TCTCCATGGTGCATCAGCAC   |     4 |
|  2 | A1BG_TGGAAGTCCACTCCACTCAG   |   487 |
|  3 | A1CF_ACAGGAAGAATTCAGTTATG   |   512 |
|  4 | A1CF_AGTTATGTTAGGTATACCCG   |   589 |
|  5 | A1CF_CTTCATTTCCCAGCCACCAA   |  1037 |
|  6 | A1CF_GAATTCAACAATATCAAACC   |   118 |
|  7 | A2ML1_AAATACTTACTGGTATCGAG  |    80 |
|  8 | A2ML1_ACCCTGGTTACTGATAACAA  |     8 |
|  9 | A2ML1_GGTGGATTATTACATCGACC  |   170 |


# ---
article_title: Genome-Scale CRISPR Screen Reveals Neddylation to Contribute to Cisplatin Resistance of Testicular Germ Cell Tumours

dataset_filename: funkeGenomescaleCRISPRScreen2023_37024667_41416_2023_2247_MOESM3_ESM_7_NGS_Cand_list_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: RefSeq ID

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['None']

sample:
|    | Gene Name   | Refseq ID    | Target sequence sgRNA   |   Read counts |   Unnamed: 4 |
|---:|:------------|:-------------|:------------------------|--------------:|-------------:|
|  0 | TRAP1       | NM_016292    | AGGGCGACGGGCCTTGCGCG    |       8044032 |          nan |
|  1 | TK1         | NM_003258    | CGTGCGTCCCTCTGTTTATA    |       7803883 |          nan |
|  2 | DDB1        | NM_001923    | GCCTTCGGATGTGGCGGATG    |       6841409 |          nan |
|  3 | ENSA        | NM_207047    | TGCTTTGGCGCTGGTTAGTT    |       1590607 |          nan |
|  4 | IFNW1       | NM_002177    | CCCAGGTACTCAGGGGGCTG    |        656648 |          nan |
|  5 | AKT3        | NM_001206729 | CAGAAGAATCGCTTGAACCT    |        531280 |          nan |
|  6 | PIR         | NM_003662    | GAGACTGAGCCACCTGTCCC    |        499632 |          nan |
|  7 | AK1         | NM_000476    | AAAGTGTAGTTCGCGTGTGT    |        496415 |          nan |
|  8 | KIF27       | NM_001271928 | CGCGTTGGTGGGACACAACT    |        489821 |          nan |
|  9 | GMNC        | NM_001146686 | GTCTTCTGCTGGGTGACTCC    |        483723 |          nan |


# ---
article_title: Genome-Wide CRISPR/Cas9 Screen Identifies Etoposide Response Modulators Associated with Clinical Outcomes in Pediatric AML

dataset_filename: nguyenGenomewideCRISPRCas92022_36111891_BLOODA_ADV-2022-007934-mmc1_Supplemental Table S0_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Etoposide', 'Vehicle control', 'Day 0', 'Day 4', 'Day 12', 'Day 18']

sample:
|    | Supplemental Table S0. Genome-wide Etoposide score (etoposide vs. vehicle controls) and essential score (vehicle control vs. day 0 replicates) at day 4 (early), day 12 (intermediate), and day 18 (late) using MAGeCK-RRA   | Unnamed: 1         | Unnamed: 2    | Unnamed: 3     | Unnamed: 4   | Unnamed: 5         | Unnamed: 6     | Unnamed: 7      | Unnamed: 8    | Unnamed: 9         | Unnamed: 10    | Unnamed: 11     | Unnamed: 12   |
|---:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------|:--------------|:---------------|:-------------|:-------------------|:---------------|:----------------|:--------------|:-------------------|:---------------|:----------------|:--------------|
|  0 | Gene                                                                                                                                                                                                                         | Day4.Etop.Score    | Day4.Etop.FDR | Day4.Ess.Score | Day4.Ess.FDR | Day12.Etop.Score   | Day12.Etop.FDR | Day12.Ess.Score | Day12.Ess.FDR | Day18.Etop.Score   | Day18.Etop.FDR | Day18.Ess.Score | Day18.Ess.FDR |
|  1 | A1BG                                                                                                                                                                                                                         | -0.982757915452354 | 0.972306      | 0.20568        | 0.853545     | 1.7519285990808    | 0.752928       | -1.7533         | 0.099903      | -1.33446591156667  | 0.933154       | -0.036103       | 1             |
|  2 | A1CF                                                                                                                                                                                                                         | 0.341254264638847  | 1             | 0.092881       | 1            | -0.483398528473466 | 1              | 0.08054         | 1             | 1.19622888552981   | 0.960285       | 0.22909         | 1             |
|  3 | A2M                                                                                                                                                                                                                          | -0.621456752838356 | 0.998435      | 0.1778         | 0.970935     | -0.715817505363098 | 1              | 0.56573         | 0.82638       | 1.07435277487943   | 0.979498       | 0.18195         | 0.97963       |
|  4 | A2ML1                                                                                                                                                                                                                        | 0.726373793727041  | 1             | -0.10209       | 0.999976     | 0.3131939780091    | 1              | 0.16833         | 1             | -0.995205889611288 | 0.978192       | 0.13114         | 1             |
|  5 | A3GALT2                                                                                                                                                                                                                      | -1.3992121814927   | 0.944459      | 0.10088        | 0.807535     | -1.35563035163587  | 0.953337       | -0.09529        | 0.891131      | -1.68468053474125  | 0.878414       | 0.25796         | 0.934135      |
|  6 | A4GALT                                                                                                                                                                                                                       | -0.118409943696597 | 1             | -0.27607       | 0.970014     | -0.185785710336726 | 1              | -0.11946        | 1             | -1.06614842205704  | 0.971275       | 0.10813         | 0.988299      |
|  7 | A4GNT                                                                                                                                                                                                                        | -0.591098917645857 | 1             | 0.099189       | 1            | 1.14783730360741   | 0.97735        | -0.26982        | 1             | 0.322557307344802  | 1              | -0.13242        | 1             |
|  8 | AAAS                                                                                                                                                                                                                         | 0.544820024454471  | 1             | -0.04374       | 0.999976     | 0.443685435664653  | 1              | 0.27547         | 1             | 1.91058948901646   | 0.77001        | -0.11169        | 1             |
|  9 | AACS                                                                                                                                                                                                                         | -0.994219100045348 | 0.970031      | -0.14916       | 0.999976     | -0.602651033823409 | 1              | 0.54512         | 0.774026      | 2.27302841631712   | 0.63101        | 0.54893         | 0.391172      |


# ---
article_title: A Whole-genome CRISPR Screen Identifies a Role of MSH2 in Cisplatin-mediated Cell Death in Muscle-invasive Bladder Cancer

dataset_filename: goodspeedWholegenomeCRISPRScreen2019_30414698_NIHMS1510205-supplement-1_Supplementary Table 6_t1

header_row: 1

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['baseMean1', 'baseMean2', 'baseMean3']

sample:
|    | Supplementary Table 6. Normalized screen data and DeSEQ2 output.    | Unnamed: 1   | Unnamed: 2   | Unnamed: 3   | Unnamed: 4   | Unnamed: 5   | Unnamed: 6   | Unnamed: 7   | Unnamed: 8   | Unnamed: 9   | Unnamed: 10   | Unnamed: 11   | Unnamed: 12   | Unnamed: 13        | Unnamed: 14           |
|---:|:--------------------------------------------------------------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:-------------|:--------------|:--------------|:--------------|:-------------------|:----------------------|
|  0 | Gene                                                                | baseMean1    | Log2_FC1     | pval1        | baseMean2    | Log2_FC2     | pval2        | baseMean3    | Log2_FC3     | pval3        | chi           | combined_p    | combined_q    | Meidan_Log2_FC     | Combined_mean_Log2_FC |
|  1 | MSH2                                                                | 957.2566405  | 5.67843972   | 1.83e-11     | 1975.317814  | 5.116799119  | 1.77e-11     | 1072.351948  | 4.822651573  | 2.38e-07     | 187           | 1.1e-37       | 1.83e-33      | 5.11679911917568   | 5.25058214863253      |
|  2 | MLH1                                                                | 1825.532851  | 3.686081025  | 7.95e-08     | 581.8766207  | 3.253182804  | 0.000102254  | 67.24306563  | 2.459707092  | 0.01856056   | 85.2          | 3e-16         | 4.97e-12      | 3.25318280400865   | 3.21759138327125      |
|  3 | FAM89B                                                              | 34.07272074  | 0.309519502  | 1            | 35.94595749  | -4.138685456 | 2.42e-06     | 418.9378725  | 3.69554062   | 5.49e-07     | 78.9          | 6.03e-15      | 9.98e-11      | 0.309519502289633  | 2.2481295828933       |
|  4 | XPC                                                                 | 242.0108639  | 3.590859153  | 0.000204781  | 510.0933263  | 2.508409082  | 0.004336734  | 3303.070108  | 3.92785658   | 1.8e-06      | 78.4          | 7.64e-15      | 1.27e-10      | 3.59085915261611   | 3.45762272302659      |
|  5 | PMS2                                                                | 185.5368054  | 1.833034834  | 0.053103685  | 38.94761474  | 0.3171524    | 0.74264162   | 2553.802775  | 4.296609819  | 2.46e-10     | 73.2          | 9e-14         | 1.49e-09      | 1.83303483405184   | 3.02743733911321      |
|  6 | LIMS2                                                               | 673.1455423  | 1.858821731  | 0.027029046  | 590.5748256  | 3.328626576  | 1.08e-05     | 16.50084367  | -3.849897954 | 0.000112886  | 69.6          | 4.94e-13      | 8.18e-09      | 1.85882173062717   | 2.19566302306935      |
|  7 | ANGPTL7                                                             | 96.50486054  | 2.140099218  | 0.017934723  | 345.7182685  | 2.322203294  | 0.01669946   | 44.92332293  | -4.145266276 | 5.96e-07     | 64.8          | 4.74e-12      | 7.85e-08      | 2.14009921801679   | 1.65769956156555      |
|  8 | PDCD5                                                               | 254.7913874  | 2.575088681  | 0.008850105  | 36.63127137  | -4.261644214 | 3.52e-07     | 277.6747389  | 1.705177402  | 0.080823638  | 63.8          | 7.58e-12      | 1.26e-07      | 1.70517740226895   | 1.62790249697668      |
|  9 | KLHL41                                                              | 43.48558655  | -4.264110805 | 7.93e-07     | 24.21294152  | -0.848017794 | 0.360969876  | 259.2117599  | 2.760069874  | 0.000859159  | 63.8          | 7.58e-12      | 1.26e-07      | -0.848017793818098 | 1.29902629084657      |


# ---
article_title: {{MND1}} and {{PSMC3IP}} Control {{PARP}} Inhibitor Sensitivity in Mitotic Cells

dataset_filename: zelceskiMND1PSMC3IPControl2023_37163373_1-s2.0-S2211124723004953-mmc2_ST1B_t1

header_row: 0

sgRNA_sequence: Own-column

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Z-score/Z-ratios']

conditions: ['olap', 'talaz']

sample:
|    | GeneID   | Ensembl gene ID   | sgRNA                                           |   CRISPRn Viability Z score (olap and talaz screen) |   CRISPRn Drug Effect Z score (olap screen) |   CRISPRn Drug Effect Z score (talaz screen) | Essentiality   |
|---:|:---------|:------------------|:------------------------------------------------|----------------------------------------------------:|--------------------------------------------:|---------------------------------------------:|:---------------|
|  0 | A1BG     | ENSG00000121410   | A1BG_CCDS12976.1_ex3_19:58862927-58862950:-_5-1 |                                            0.668991 |                                    0.705699 |                                    -0.505318 | non-essential  |
|  1 | A1BG     | ENSG00000121410   | A1BG_CCDS12976.1_ex5_19:58864367-58864390:-_5-5 |                                           -1.15848  |                                    0.467419 |                                     0.022966 | non-essential  |
|  2 | A1BG     | ENSG00000121410   | A1BG_CCDS12976.1_ex4_19:58863655-58863678:+_5-2 |                                           -0.090289 |                                    0.260061 |                                     0.814947 | non-essential  |
|  3 | A1BG     | ENSG00000121410   | A1BG_CCDS12976.1_ex4_19:58863697-58863720:-_5-3 |                                           -1.12321  |                                    0.215305 |                                     0.754094 | non-essential  |
|  4 | A1BG     | ENSG00000121410   | A1BG_CCDS12976.1_ex4_19:58863866-58863889:+_5-4 |                                           -2.07258  |                                   -0.373944 |                                     0.446271 | non-essential  |
|  5 | A1CF     | ENSG00000148584   | A1CF_CCDS7241.1_ex9_10:52603844-52603867:-_5-5  |                                           -1.01605  |                                   -0.231366 |                                    -0.244258 | non-essential  |
|  6 | A1CF     | ENSG00000148584   | A1CF_CCDS7241.1_ex7_10:52596023-52596046:+_5-3  |                                            0.135322 |                                    0.266967 |                                     1.0657   | non-essential  |
|  7 | A1CF     | ENSG00000148584   | A1CF_CCDS7241.1_ex6_10:52588014-52588037:-_5-1  |                                            0.264492 |                                   -0.483242 |                                    -0.449798 | non-essential  |
|  8 | A1CF     | ENSG00000148584   | A1CF_CCDS7241.1_ex9_10:52603761-52603784:+_5-4  |                                            2.42248  |                                    2.45976  |                                     2.645    | non-essential  |
|  9 | A1CF     | ENSG00000148584   | A1CF_CCDS7241.1_ex7_10:52595962-52595985:-_5-2  |                                            4.81517  |                                   -4.77792  |                                     3.44597  | non-essential  |


# ---
article_title: Genome-Wide CRISPR Screen Identifies ELP5 as a Determinant of Gemcitabine Sensitivity in Gallbladder Cancer

dataset_filename: xuGenomewideCRISPRScreen2019_31792210_41467_2019_13420_MOESM4_ESM_reads out counts_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['Survive', 'Input']

sample:
|    | GeneID   |   Survive_1 |   Survive_2 |   Input_1 |   Input_2 |
|---:|:---------|------------:|------------:|----------:|----------:|
|  0 | A1BG     |           0 |           0 |         0 |         0 |
|  1 | A1BG     |         396 |         219 |       153 |       154 |
|  2 | A1BG     |         286 |         362 |       226 |       262 |
|  3 | A1BG     |         325 |         237 |       275 |       266 |
|  4 | A1CF     |         128 |         116 |       129 |       132 |
|  5 | A1CF     |          95 |          66 |       198 |       180 |
|  6 | A1CF     |         226 |         214 |       234 |       215 |
|  7 | A1CF     |         437 |         257 |       261 |       238 |
|  8 | A2M      |         104 |         125 |        87 |        99 |
|  9 | A2M      |         222 |         231 |       171 |       120 |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_Maphosphamide_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Maphosphamide']

sample:
|    | id     |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:-------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | TP53   |     6 |     1       |       1       |         1 |      20466 |               0 |   2.1621  |  1.7926e-19 |    2.4189e-07 |   0.00033 |          1 |               6 |   2.1621  |
|  1 | HEATR3 |     6 |     1       |       1       |         1 |      20460 |               0 |   1.2388  |  3.7987e-12 |    2.4189e-07 |   0.00033 |          2 |               6 |   1.2388  |
|  2 | KAT5   |     6 |     0.57038 |       0.71191 |         1 |      14108 |               1 |   1.2079  |  5.6488e-12 |    2.4189e-07 |   0.00033 |          3 |               5 |   1.2079  |
|  3 | HNRNPM |     6 |     0.90892 |       0.91167 |         1 |      18329 |               1 |   1.0584  |  8.3701e-11 |    2.4189e-07 |   0.00033 |          4 |               5 |   1.0584  |
|  4 | PMAIP1 |     6 |     1       |       1       |         1 |      20465 |               0 |   0.84946 |  8.5806e-11 |    2.4189e-07 |   0.00033 |          5 |               6 |   0.84946 |
|  5 | ACSL4  |     6 |     0.98696 |       0.98692 |         1 |      20020 |               0 |   1.4665  |  1.2569e-10 |    2.4189e-07 |   0.00033 |          6 |               5 |   1.4665  |
|  6 | EPC2   |     6 |     1       |       1       |         1 |      20463 |               0 |   1.058   |  1.7193e-10 |    2.4189e-07 |   0.00033 |          7 |               6 |   1.058   |
|  7 | FLCN   |     6 |     0.99998 |       0.99998 |         1 |      20447 |               0 |   1.2395  |  2.0997e-10 |    2.4189e-07 |   0.00033 |          8 |               6 |   1.2395  |
|  8 | RPL11  |     6 |     0.98995 |       0.98995 |         1 |      20106 |               0 |   1.5546  |  3.5135e-10 |    2.4189e-07 |   0.00033 |          9 |               5 |   1.5546  |
|  9 | EBF1   |     6 |     0.17108 |       0.3484  |         1 |       6924 |               1 |   0.84909 |  7.2146e-10 |    2.4189e-07 |   0.00033 |         10 |               3 |   0.84909 |


# ---
article_title: A Genome-Scale CRISPR/Cas9 Knockout Screening Reveals SH3D21 as a Sensitizer for Gemcitabine

dataset_filename: masoudiGenomescaleCRISPRCas92019_31844142_41598_2019_55893_MOESM3_ESM_Sheet1_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['Gemcitabine']

sample:
|    | Gene     |   # Hairpins |       NES |   Gene rank |   p-value |   p-value rank |
|---:|:---------|-------------:|----------:|------------:|----------:|---------------:|
|  0 | KIAA0196 |            6 | 0.0005764 |           1 |    0.0001 |              4 |
|  1 | SLFN14   |            6 | 0.004255  |           2 |    0.0001 |              2 |
|  2 | MTRNR2L2 |            6 | 0.004413  |           3 |    0.0001 |              1 |
|  3 | DGCR6    |            6 | 0.005611  |           4 |    0.0004 |              6 |
|  4 | TMSB4Y   |            4 | 0.005614  |           5 |    0.0001 |              3 |
|  5 | OR2A25   |            6 | 0.006507  |           6 |    0.0004 |              7 |
|  6 | DSTN     |            6 | 0.006534  |           7 |    0.0004 |              8 |
|  7 | RPTN     |            6 | 0.007925  |           8 |    0.0005 |             12 |
|  8 | ZNF446   |            6 | 0.008105  |           9 |    0.0005 |             10 |
|  9 | RPL3     |            6 | 0.008929  |          10 |    0.0005 |             11 |


# ---
article_title: {{CCNE1}} Amplification Is Synthetic Lethal with {{PKMYT1}} Kinase Inhibition

dataset_filename: galloCCNE1AmplificationSynthetic2022_35444283_41586_2022_4638_MOESM3_ESM_ReadCounts_TKOv3 CCNE1 SL_t1

header_row: 0

sgRNA_sequence: Concatenated

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per sgRNA

metrics: ['Raw count']

conditions: ['TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0', 'TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18A', 'TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18B', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T0', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18A', 'TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18B']

sample:
|    | Read counts for screens using TKOv3   | Unnamed: 1   | Unnamed: 2                                     | Unnamed: 3                                       | Unnamed: 4                                       | Unnamed: 5                        | Unnamed: 6                          | Unnamed: 7                          |
|---:|:--------------------------------------|:-------------|:-----------------------------------------------|:-------------------------------------------------|:-------------------------------------------------|:----------------------------------|:------------------------------------|:------------------------------------|
|  0 | sgRNA                                 | Gene         | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T0 | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18A | TKOv3_RPE1.hTERT.Cas9.dTP53_CCNE1oe_Clone#2_T18B | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T0 | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18A | TKOv3_RPE1.hTERT.Cas9.dTP53_WT_T18B |
|  1 | TKOv3_A1BG_exon1_4                    | A1BG         | 1501                                           | 1029                                             | 240                                              | 451                               | 366                                 | 440                                 |
|  2 | TKOv3_A1BG_exon3_3                    | A1BG         | 1404                                           | 703                                              | 463                                              | 404                               | 352                                 | 404                                 |
|  3 | TKOv3_A1BG_exon4_2                    | A1BG         | 535                                            | 401                                              | 266                                              | 146                               | 66                                  | 179                                 |
|  4 | TKOv3_A1BG_exon5_1                    | A1BG         | 708                                            | 400                                              | 219                                              | 236                               | 111                                 | 252                                 |
|  5 | TKOv3_A1CF_exon1_4                    | A1CF         | 366                                            | 281                                              | 288                                              | 166                               | 141                                 | 95                                  |
|  6 | TKOv3_A1CF_exon2_3                    | A1CF         | 398                                            | 395                                              | 285                                              | 154                               | 91                                  | 121                                 |
|  7 | TKOv3_A1CF_exon4_2                    | A1CF         | 720                                            | 1301                                             | 818                                              | 206                               | 235                                 | 196                                 |
|  8 | TKOv3_A1CF_exon7_1                    | A1CF         | 1471                                           | 1336                                             | 919                                              | 450                               | 403                                 | 378                                 |
|  9 | TKOv3_A2ML1_exon1_1                   | A2ML1        | 845                                            | 974                                              | 878                                              | 340                               | 263                                 | 202                                 |


# ---
article_title: Mutational and Functional Genetics Mapping of Chemotherapy Resistance Mechanisms in Relapsed Acute Lymphoblastic Leukemia

dataset_filename: oshimaMutationalFunctionalGenetics2020_33796864_NIHMS1673078-supplement-Supplementary_Table_7_AraC_t1

header_row: 0

sgRNA_sequence: None

gene_identifier: HGNC gene symbol

dataset_type: Data

statistic_aggregation: Per gene

metrics: ['MAGeCK statistics']

conditions: ['AraC']

sample:
|    | id      |   num |   neg_score |   neg_p-value |   neg_fdr |   neg_rank |   neg_goodsgrna |   neg_lfc |   pos_score |   pos_p-value |   pos_fdr |   pos_rank |   pos_goodsgrna |   pos_lfc |
|---:|:--------|------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|------------:|--------------:|----------:|-----------:|----------------:|----------:|
|  0 | DCK     |     6 |           1 |             1 |         1 |      20466 |               0 |    3.3873 |  1.1002e-22 |    2.4189e-07 |  0.000177 |          1 |               6 |    3.3873 |
|  1 | SLC29A1 |     6 |           1 |             1 |         1 |      20462 |               0 |    3.2527 |  4.7006e-15 |    2.4189e-07 |  0.000177 |          2 |               6 |    3.2527 |
|  2 | NOP16   |     6 |           1 |             1 |         1 |      20465 |               0 |    1.5237 |  8.6231e-13 |    2.4189e-07 |  0.000177 |          3 |               6 |    1.5237 |
|  3 | IARS2   |     6 |           1 |             1 |         1 |      20464 |               0 |    1.7807 |  2.1605e-12 |    2.4189e-07 |  0.000177 |          4 |               6 |    1.7807 |
|  4 | FARSA   |     6 |           1 |             1 |         1 |      20463 |               0 |    1.2936 |  1.5707e-11 |    2.4189e-07 |  0.000177 |          5 |               6 |    1.2936 |
|  5 | EIF2S3  |     6 |           1 |             1 |         1 |      20461 |               0 |    1.3894 |  5.0349e-11 |    2.4189e-07 |  0.000177 |          6 |               6 |    1.3894 |
|  6 | HSPE1   |     5 |           1 |             1 |         1 |      20460 |               0 |    1.3416 |  7.9797e-11 |    2.4189e-07 |  0.000177 |          7 |               5 |    1.3416 |
|  7 | TRMT112 |     6 |           1 |             1 |         1 |      20459 |               0 |    1.8207 |  8.008e-11  |    2.4189e-07 |  0.000177 |          8 |               6 |    1.8207 |
|  8 | PELO    |     6 |           1 |             1 |         1 |      20458 |               0 |    1.2798 |  8.7076e-11 |    2.4189e-07 |  0.000177 |          9 |               6 |    1.2798 |
|  9 | IMPDH2  |     6 |           1 |             1 |         1 |      20457 |               0 |    1.5845 |  9.6066e-11 |    2.4189e-07 |  0.000177 |         10 |               6 |    1.5845 |


