library(tidyverse)

df <- readr::read_tsv(
  snakemake@input[[1]],
  col_names = c("chrom", "start", "end", "gene", "mapq", "score"),
  col_types = "ciiccc"
)

df_expanded <- df %>%
  mutate(gene = str_split(gene, ","),
         mapq = str_split(mapq, ","),
         score = str_split(score, ",")) %>%
  unnest(c(gene, mapq, score)) %>%
  mutate(mapq = as.numeric(mapq),
         score = as.numeric(str_replace(score, "AS:i:", "")),
         gene = str_extract(gene, ".*\\|(.*)\\*[0-9]+\\|.*", 1),
         .locus = str_sub(gene, 1, 3),
         .chrom = str_extract(chrom, "(chr[0-9]+)_?.*", 1),
         scorefrac = score / (end - start)) %>%
  mutate(correct_map = case_when(.locus == "TRA" & .chrom == "chr14" ~ TRUE,
                                 .locus == "TRD" & .chrom == "chr14" ~ TRUE,
                                 .locus == "TRB" & .chrom == "chr7" ~ TRUE,
                                 .locus == "TRG" & .chrom == "chr7" ~ TRUE,
                                 .locus == "IGK" & .chrom == "chr2" ~ TRUE,
                                 .locus == "IGL" & .chrom == "chr22" ~ TRUE,
                                 .locus == "IGH" & .chrom == "chr14" & start > 90000000 ~ TRUE,
                                 TRUE ~ FALSE)) %>%
  select(-.locus, -.chrom)

df_totally_wrong <- df_expanded %>%
  filter(!correct_map)

df_correct <- df_expanded %>%
  filter(correct_map)

df_allele <- df_correct %>%
  group_by(chrom, start, end, gene) %>%
  summarize(n_alleles = n(),
            scoreS = sd(score),
            score = mean(score),
            scorefracS = sd(scorefrac),
            scorefrac = mean(scorefrac),
            mapqS = sd(mapq),
            mapq = mean(mapq),
            .groups = "drop")

# eliminate locations for which the highest match score is 20%
# and all genes whose highest mapq is 10
df_max <- df_allele %>%
  group_by(chrom, start, end) %>%
  mutate(scoremax_location = max(scorefrac),
         mapqmax_location = max(mapq)) %>%
  group_by(gene) %>%
  mutate(mapqmax_gene = max(mapq),
         scoremax_gene = max(scorefrac)) %>%
  ungroup()

## df_max %>%
##   filter(scoremax_location < 0.2)

## df_max %>%
##   filter(mapqmax_location < 5)

## df_max %>%
##   filter(scoremax_gene < 0.2)

## df_max %>%
##   filter(mapqmax_gene < 5)

# find genes that map to only one location, which in turn only has one gene
df_counted <- df_max %>%
  filter(scoremax_location > 0.2 &
           mapqmax_location > 5 &
           scoremax_gene > 0.2 &
           mapqmax_gene > 5) %>%
  select(-matches("max")) %>%
  group_by(chrom, start, end) %>%
  arrange(score) %>%
  mutate(n_genes = n()) %>%
  ungroup() %>%
  group_by(gene) %>%
  arrange(mapq) %>%
  mutate(n_locations = n()) %>%
  ungroup()

df_single_unique <- df_counted %>%
  filter(n_locations == 1 & n_genes == 1)

# find genes that rank the best in terms of alignment score and map quality at
# a given location and within a given gene; use this to eliminate gene/location
# choices
df_diff <- df_counted %>%
  anti_join(df_single_unique, by = c("chrom", "start", "end")) %>%
  anti_join(df_single_unique, by = "gene") %>%
  group_by(chrom, start, end) %>%
  arrange(score) %>%
  mutate(scorediff_loc = scorefrac - lag(scorefrac),
         score_loc_rank = n() - row_number()) %>%
  arrange(mapq) %>%
  mutate(mapqdiff_loc = mapq - lag(mapq),
         mapq_loc_rank = n() - row_number()) %>%
  group_by(gene) %>%
  arrange(mapq) %>%
  mutate(mapqdiff_gene = mapq - lag(mapq),
         mapq_gene_rank = n() - row_number()) %>%
  arrange(scorefrac) %>%
  mutate(scorediff_gene = scorefrac - lag(scorefrac),
         score_gene_rank = n() - row_number()) %>%
  ungroup() %>%
  replace_na(list(scorediff_loc = 0,
                  scorediff_gene = 0,
                  mapqdiff_loc = 0,
                  mapqdiff_gene = 0)) %>%
  relocate(scorefrac, mapq, matches("diff"), matches("rank"), .after = gene)

df_high_map_score <- df_diff %>%
  filter(scorediff_loc > 0.25 & score_loc_rank == 0 &
           scorediff_gene > 0.25 & score_gene_rank == 0 &
           mapqdiff_loc > 20 & mapq_loc_rank == 0 &
           mapqdiff_gene > 20 & mapq_gene_rank == 0) %>%
  select(-matches("n_")) %>%
  group_by(chrom, start, end) %>%
  mutate(n_genes = n()) %>%
  group_by(gene) %>%
  mutate(n_locations = n()) %>%
  ungroup()

df_high_rank_diff_unique <- df_high_map_score %>%
  filter(n_locations == 1 & n_genes == 1)

df_high_rank_unique <- df_diff %>%
  anti_join(df_high_rank_diff_unique, by = c("chrom", "start", "end")) %>%
  anti_join(df_high_rank_diff_unique, by = "gene") %>%
  filter(score_loc_rank == 0 &
           score_gene_rank == 0 &
           mapq_loc_rank == 0 &
           mapq_gene_rank == 0) %>%
  select(-matches("n_")) %>%
  group_by(chrom, start, end) %>%
  mutate(n_genes = n()) %>%
  group_by(gene) %>%
  mutate(n_locations = n()) %>%
  ungroup() %>%
  filter(n_locations == 1 & n_genes == 1)

df_unique <- df_high_rank_unique %>%
  bind_rows(df_high_rank_diff_unique) %>%
  bind_rows(df_single_unique)

df_unique %>%
  select(chrom, start, end, gene) %>%
  readr::write_tsv(snakemake@output[[1]], col_names = F)
