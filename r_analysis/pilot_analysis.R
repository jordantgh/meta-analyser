box::use(
  readxl[rxl = read_excel],
  glue[g = glue],
  dplyr[...],
  ggplot2[...]
)

df <- rxl(
  g(
    "radchemoscreens/downloads/",
    "awahRibosomalProteinS112020_32528131",
    "_NIHMS1599130-supplement-1599130_SuppList_Data_6.xlsx"
  ),
  range = "B1:D77442"
) %>% cbind(
  rxl(
    g(
      "radchemoscreens/downloads/",
      "awahRibosomalProteinS112020_32528131",
      "_NIHMS1599130-supplement-1599130_SuppList_Data_6.xlsx"
    ),
    range = "G1:G77442"
  ),
  .
)

df_counts <- df %>%
  filter(!if_any(everything(), \(x) x == 0)) %>%
  mutate_at(vars(-c("genes")), log2)

# plot the eto counts against the dmso counts
plot <- ggplot(df_counts, aes(x = DMSO, y = !!sym("Eto R1"))) +
  geom_point() +
  ylim(0, 15) +
  xlim(0, 15)

# convert counts to fraction of column mean
df <- df %>%
  mutate_at(vars(-c("genes")), \(x) x / mean(x)) %>%
  filter(DMSO != 0)


# calculate average change relative to dmso
df <- df %>%
  mutate(
    eto_avg = (!!sym("Eto R1") + !!sym("Eto R2")) / 2,
    rel_change = eto_avg / DMSO
  ) %>%
  group_by(genes) %>%
  summarise(
    dmso = mean(DMSO),
    eto_avg = mean(eto_avg),
    rel_change = mean(rel_change)
  )

# load isgs csv
isgs <- read.csv("wilson_isg_results.csv")

# filter on species column
isgs <- isgs %>%
  select(-c(Orthologous.Cluster.ID)) %>%
  filter(Species == "Homo sapiens" & Expression == "up_regulated")

# plot eto_avg against dmso
library(ggplot2)

plot <- ggplot(df, aes(x = dmso, y = eto_avg)) +
  geom_point()

# Join the df and isgs dataframes on the gene column
merged_df <- df %>% inner_join(isgs, by = c("genes" = "Gene"))

# Calculate the average rel_change for the genes present in the isgs dataframe
overall_avg_change <- mean(df$rel_change)
isg_avg_change <- mean(merged_df$rel_change)
