# Load required libraries
library(DBI)

# Connect to SQLite database
con <- dbConnect(
  RSQLite::SQLite(),
  "post_pruning_tables-70217c49-e03c-4985-8a35-84e739b8a7e0.db"
)

# Get list of all table names in the database
table_names <- dbListTables(con)

# Open Markdown file for writing
markdown_file_path <- "table_previews.md"
writeLines("", markdown_file_path)
filecon <- file(markdown_file_path, "a")

# Loop through each table
for (table_name in table_names) {
  # Write table name as H2 header in Markdown
  writeLines(paste0("## ", table_name), filecon)

  # Fetch all rows from the table (this will load the entire table into memory)
  data <- dbReadTable(con, table_name)

  # Convert data to Markdown table using knitr::kable and write to file
  writeLines(
    knitr::kable(head(data, 10), format = "markdown"),
    filecon
  )

  # Add a new line for separation
  writeLines("\n", filecon)
}

# Close Markdown file
close(filecon)

# Disconnect from the database
dbDisconnect(con)


{
  "header_row_index": "0",
  "sgRNA_sequence": "none",
  "gene_identifier": "HGNC gene symbol",
  "measurements_per_gene": "one",
  "all_applicable_metrics": [
    "p-value",
    "fdr",
    "fold-change/log fold-change",
    "rank"
  ]
}