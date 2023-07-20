library(tidyverse)

order_df <- readr::read_tsv(
  snakemake@input$order,
  col_types = "ci"
) %>%
  # dummy column
  mutate(pos_real = pos)

# TODO not DRY at alllllllllllll
read_bed <- function(path) {
  readr::read_tsv(
    path,
    col_types = "ciic",
    col_names = c("chrom", "start", "end", "gene")
  ) %>%
    mutate(gene = str_replace_all(gene, "[()]", "")) %>%
    # probably don't need orphons since most (all?) of these are outside the main
    # loci anyways
    filter(!str_detect(gene, "OR")) %>%
    arrange(chrom, start) %>%
    mutate(mapped_pos = row_number()) %>%
    full_join(order_df, by = "gene") %>%
    mutate(locus = str_sub(gene, 1, 3),
           locus = if_else(locus == "TRD", "TRA", locus),
           major = str_sub(gene, 4, 4),
           minor = str_sub(gene, 5)) %>%
    mutate(is_constant = !major %in% c("V", "D", "J", "DEL"),
           minor = if_else(is_constant, paste0(major, minor), minor),
           major = if_else(is_constant, "const", major))
}

read_bed(snakemake@input$bed) %>%
  # dummy columns
  mutate(correct_map = TRUE, hap = "unit") %>%
  readr::write_tsv(snakemake@output[[1]])
