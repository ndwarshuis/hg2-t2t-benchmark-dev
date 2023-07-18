library(tidyverse)

split_attrs <- function(s) {
  res <- str_split_1(s, ";") %>%
    map(~ as.list(str_split_1(.x, "="))) %>%
    transpose()
  tibble(key = unlist(res[[1]]), value = unlist(res[[2]]))
}

read_mapped_data <- function(order_df, mapped_path, unmapped_path, hap) {
  gff <- read_tsv(
    mapped_path,
    comment = "#",
    col_names = c("chrom", NA, NA,  "start", "end", NA, NA, NA, "attr"),
    col_types = cols(chrom = "c", start = "i", end = "i", "attr" = "c", .default = "-")
  ) %>%
    mutate(attr = map(attr, split_attrs),
           chrom = str_extract(chrom, "(chr.*?)_[^_]+", 1)) %>%
    unnest(attr) %>%
    filter(key %in% c("Name", "description", "coverage", "sequence_ID")) %>%
    pivot_wider(id_cols = c(chrom, start, end), names_from = "key", values_from = "value") %>%
    mutate(coverage = as.numeric(coverage),
           sequence_ID = as.numeric(sequence_ID)) %>%
    rename(gene = Name)

  unmapped <- read_tsv(
    unmapped_path,
    col_types = "c",
    col_names = c("gene")
  ) %>%
    mutate(gene = str_replace(gene, "gene-", ""))

  bind_rows(unmapped, gff) %>%
    mutate(locus = str_sub(gene, 1, 3)) %>%
    mutate(locus = if_else(locus == "TRD", "TRA", locus)) %>%
    mutate(correct_map = case_when(gene == "IGKDEL" & chrom == "chr2" ~ TRUE,
				   locus == "TRA" & chrom == "chr14" ~ TRUE,
                                   locus == "TRB" & chrom == "chr7" ~ TRUE,
                                   locus == "TRG" & chrom == "chr7" ~ TRUE,
                                   locus == "IGK" & chrom == "chr2" ~ TRUE,
                                   locus == "IGL" & chrom == "chr22" ~ TRUE,
                                   locus == "IGH" & chrom == "chr14" & start > 90000000 ~ TRUE,
                                   TRUE ~ FALSE),
           is_global = gene %in% c("IGK", "IGH", "IGL", "TRA", "TRB", "TRG")) %>%
    group_by(locus, correct_map, is_global) %>%
    arrange(chrom, start) %>%
    mutate(mapped_pos = row_number()) %>%
    ungroup() %>%
    select(-locus) %>%
    mutate(hap = hap)
}

order_df <- read_tsv(
  snakemake@input$order,
  col_types = "ci",
  col_names = c("gene", "pos")
)

pat_query <- read_mapped_data(
  order_df,
  snakemake@input$mapped_pat,
  snakemake@input$unmapped_pat,
  "pat"
)

mat_query <- read_mapped_data(
  order_df,
  snakemake@input$mapped_mat,
  snakemake@input$unmapped_mat,
  "mat"
)

query <- bind_rows(pat_query, mat_query)

query_global <- query %>%
  filter(is_global) %>%
  filter(correct_map) %>%
  select(gene, start, end, hap) %>%
  rename(global_start = start,
         global_end = end)

df <- query %>%
  filter(!is_global) %>%
  select(-is_global) %>%
  full_join(order_df, by = "gene") %>%
  mutate(locus = str_sub(gene, 1, 3),
	 locus = if_else(locus == "TRD", "TRA", locus),
         major = if_else(gene == "IGKDEL", "DEL", str_sub(gene, 4, 4)),
         minor = if_else(gene == "IGKDEL", "DEL", str_sub(gene, 5))) %>%
  mutate(is_constant = !major %in% c("V", "D", "J", "DEL"),
         minor = if_else(is_constant, paste0(major, minor), minor),
         major = if_else(is_constant, "const", major)) %>%
  select(-gene, -is_constant) %>%
  arrange(locus, pos) %>%
  group_by(locus) %>%
  mutate(pos_real = row_number()) %>%
  ungroup()

readr::write_tsv(query_global, snakemake@output$loci)
readr::write_tsv(df, snakemake@output$genes)
