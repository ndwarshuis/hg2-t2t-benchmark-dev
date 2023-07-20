library(tidyverse)

split_attrs <- function(s) {
  res <- str_split_1(s, ";") %>%
    map(~ as.list(str_split_1(.x, "="))) %>%
    transpose()
  tibble(key = unlist(res[[1]]), value = unlist(res[[2]]))
}

# TODO dry as the amazon
read_mapped_data <- function(mapped_path) {
  read_tsv(
    mapped_path,
    comment = "#",
    col_names = c("chrom", NA, NA,  "start", "end", NA, NA, NA, "attr"),
    col_types = cols(chrom = "c", start = "i", end = "i", "attr" = "c", .default = "-")
  ) %>%
    mutate(attr = map(attr, split_attrs),
           chrom = str_extract(chrom, "(chr.*?)_[^_]+", 1)) %>%
    unnest(attr) %>%
    filter(key %in% c("Name", "description")) %>%
    pivot_wider(id_cols = c(chrom, start, end), names_from = "key", values_from = "value") %>%
    rename(gene = Name) %>%
    mutate(is_global = gene %in% c("IGK", "IGH", "IGL", "TRA", "TRB", "TRG")) %>%
    arrange(chrom, start)
}

query <- read_mapped_data(snakemake@input[[1]])

query_global <- query %>%
  filter(is_global) %>%
  select(gene, start, end) %>%
  # dummy column for later
  mutate(hap = "unit") %>%
  rename(global_start = start,
         global_end = end)

df <- query %>%
  filter(!is_global) %>%
  select(-is_global) %>%
  mutate(locus = str_sub(gene, 1, 3),
         locus = if_else(locus == "TRD", "TRA", locus),
         major = if_else(gene == "IGKDEL", "DEL", str_sub(gene, 4, 4)),
         minor = if_else(gene == "IGKDEL", "DEL", str_sub(gene, 5))) %>%
  mutate(is_constant = !major %in% c("V", "D", "J", "DEL"),
         minor = if_else(is_constant, paste0(major, minor), minor),
         major = if_else(is_constant, "const", major)) %>%
  select(-gene, -is_constant)

readr::write_tsv(query_global, snakemake@output$loci)
readr::write_tsv(df, snakemake@output$genes)
