library(tidyverse)

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  mutate(major = str_sub(gene, 4, 4)) %>%
  filter(major %in% c("V", "D", "J")) %>%
  filter(gene != "IGHD") %>%
  select(-major)

read_cov <- function(path) {
  read_tsv(path, col_types = "ciiiidddd") %>%
    mutate(id = basename(dirname(path)))
}

snakemake@input %>%
  map(read_cov) %>%
  list_rbind() %>%
  mutate(gene = str_extract(`#rname`, ".+\\|(.*)\\|.+", 1)) %>%
  select(-`#rname`) %>%
  separate(col = gene, into = c("gene", "allele"), sep = "\\*") %>%
  mutate(gene = str_replace_all(gene, "[()]", "")) %>%
  left_join(order_df, by = "gene") %>%
  mutate(locus = str_sub(gene, 1, 3),
         major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5),
         locus = if_else(locus == "TRD", "TRA", locus)) %>%
  filter(numreads > 0) %>%
  filter(!is.na(pos)) %>%
  filter(coverage > 50) %>%
  # for the few remaining rows that have two alleles, keep the highest mapQ
  group_by(id, gene) %>%
  filter(meanmapq == max(meanmapq)) %>%
  filter(row_number() == 1) %>%
  ungroup()

cov_df %>%
  ggplot(aes(pos, meandepth, fill = major)) +
  geom_col() +
  facet_wrap(c("id", "locus"), ncol = 6, scale = "free_x")
ggsave(snakemake@output[[1]])

## gene_df <- cov_df %>%
##   pivot_wider(id_cols = id,
##               names_from = gene,
##               values_from = meandepth,
##               values_fill = 0) %>%
##   select(id, matches("IGKJ"))

## res <- gene_df %>%
##   select(-id) %>%
##   prcomp(scale = TRUE)

## gene_df %>%
##   select(-id) %>%
##   as.matrix() %>%
##   {. %*% res$rotation} %>%
##   as_tibble() %>%
##   mutate(id = gene_df$id) %>%
##   ggplot(aes(PC1, PC2)) +
##   geom_point() +
##   geom_label(aes(label = id))
