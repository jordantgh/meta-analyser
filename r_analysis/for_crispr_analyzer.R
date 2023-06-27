box::use(
  readxl[rxl = read_excel],
  glue[g = glue],
  dplyr[...],
  tidyr[...],
  edgeR[...],
  ggplot2[...]
)

downloads <- g("{here::here()}/radchemoscreens/downloads")
article <- "awahRibosomalProteinS112020_32528131"
file <- "NIHMS1599130-supplement-1599130_SuppList_Data_6.xlsx"
full_path <- g("{downloads}/{article}_{file}")

df <- rxl(full_path, range = "A1:D77442") %>%
  rename(
    eto_r1 = `Eto R1`,
    eto_r2 = `Eto R2`
  )

# read the first data frame
df1 <- rxl(full_path, range = "A1:D77442")

# read the second data frame
df2 <- rxl(full_path, range = "G1:I77442") %>%
  mutate(genes = stringr::str_extract(sgRNA, "^[^_]*"))

# add an additional column that indicates the order of gRNA for each gene
df1 <- df1 %>%
  group_by(genes) %>%
  mutate(gRNA_order = row_number())
df2 <- df2 %>%
  group_by(genes) %>%
  mutate(gRNA_order = row_number())

# join the two data frames
df3 <- df1 %>%
  left_join(df2, by = c("genes", "gRNA_order")) %>%
  ungroup() %>%
  select(-c(genes, `Puromycin R1`, `Puromycin R2`, gRNA_order))

head(df3)

