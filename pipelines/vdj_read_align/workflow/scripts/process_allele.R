suppressMessages(library(tidyverse))

get_merge_end <- function(last, x) {
  if (x[[1]] > last[[2]]) {
    x
  } else {
    last
  }
}

map_params <- snakemake@params[["map_params"]]

optimize_order <- function(df, gap_weight, align_weight) {
  xs <- df %>%
    group_by(pos_index) %>%
    group_map(~ list(.y$pos_index[[1]], .x)) %>%
    set_names(map(., ~ .x[[1]])) %>%
    map(~ .x[[2]])
  branches <- accumulate(xs, ~ .x * nrow(.y), .init = 1)[-1] %>%
    as.vector()
  ngenes <- length(xs)
  ncombos <- last(branches)
  get_combos <- function(k) {
    xs %>%
        map(~ as.vector(.x[[k]])) %>%
        map2(branches, ~ rep(.x, .y / length(.x), each = ncombos / .y)) %>%
        flatten_dbl() %>%
        matrix(nrow = ncombos)
  }
  paths <- get_combos("gene_pos")
  dups <- apply(paths, 1, function(x) any(duplicated(x)))
  ncombos_nodup <- sum(!dups)
  unique_paths <- paths[!dups, ] %>%
    matrix(nrow = ncombos_nodup)
  aligns <- get_combos("scorefrac")[!dups, ] %>%
    matrix(nrow = ncombos_nodup)
  genes <- get_combos("gene_index")[!dups, ] %>%
    matrix(nrow = ncombos_nodup)
  gaps <- (unique_paths[, -ncol(unique_paths)] - unique_paths[, -1]) %>%
    abs() %>%
    magrittr::add(-1) %>%
    matrix(nrow = ncombos_nodup)
  gap_scores <- apply(gaps, 1, sum)
  align_scores <- apply(1 - aligns, 1, sum)
  scores <- gap_scores * gap_weight + align_scores * align_weight
  best_combo <- detect_index(scores == min(scores), ~ .x == TRUE)
  tibble(gene_index = genes[best_combo, ]) %>%
    mutate(pos_index = as.integer(names(xs)))
}

## BLACKLIST <- c("IGHV(II)-30-21-RSS*01")

CONSTANTS <- c(
  paste0("IGHC", c(paste0("G", c(1:4, "P")), "E", "M", "D", "A1", "A2", "EP1", "EP2")),
  "IGKC",
  "IGKCDEL",
  paste0("IGLC", 1:7),
  "TRAC",
  "TRDC",
  "TRBC1",
  "TRBC2",
  "TRGC1",
  "TRGC2"
) %>% tibble(
  gene = .,
  allele = "01",
  functionality = "",
  gene_orientation = "+"
)

order_df <- read_tsv(
  snakemake@input[["order"]],
  col_types = "ci"
) %>%
  rename(gene_pos = pos)

dups_df <- read_tsv(
  snakemake@input[["dups"]],
  col_types = "icc"
) %>%
  select(new, old) %>%
  rename(gene = old, new_gene = new)

headers_df <- read_tsv(
  snakemake@input[["headers"]],
  col_names = c("gene", "allele", "functionality", "gene_orientation"),
  col_types = "c"
) %>%
  bind_rows(CONSTANTS) %>%
  inner_join(order_df, by = "gene") %>%
  filter(functionality != "ORF") %>%
  select(-functionality) %>%
  mutate(gene_revcomp = gene_orientation == "-") %>%
  select(-gene_orientation)

reads_df <- read_tsv(
  snakemake@input[["reads"]],
  col_names = FALSE,
  col_types = c(X1 = "c", X2 = "i", .default = "-")
) %>%
  rename(rname = X1, rlen = X2) %>%
  mutate(rname = paste0("@", rname))

tags_df <- read_lines(snakemake@input[["tags"]]) %>%
  tibble(tags = .) %>%
  mutate(id = row_number(),
         tags = str_split(tags, "\t")) %>%
  unnest(tags) %>%
  mutate(key = str_extract(tags, "^(.*):.*:.*", 1),
         value = str_extract(tags, ".*:.*:(.*)$", 1)) %>%
  filter(key %in% c("AS", "NM")) %>%
  mutate(value = as.numeric(value)) %>%
  pivot_wider(id_cols = id, names_from = key, values_from = value) %>%
  rename(edit_distance = NM,
         score = AS) %>%
  select(-id)

sam_df <- readr::read_tsv(
  snakemake@input[["sam"]],
  col_names = c(
    "gene",
    "flag",
    "rname",
    "pos",
    "mapq",
    "cigar",
    "rnext",
    "pnext",
    "tlen",
    "seq",
    "qual"
  ),
  col_types = "cicidcciicc"
) %>%
  bind_cols(tags_df)

df <- sam_df %>%
  left_join(reads_df, by = "rname") %>%
  mutate(is_reversed = (flag %% 32) %/% 16 == 1,
         gene_index = as.integer(factor(gene)),
         seqlen = str_length(seq),
         end = map2_dbl(pos + seqlen, rlen, min),
         scorefrac = score / seqlen,
         left_clip = as.integer(str_extract(cigar, "^([0-9]+)S.*", 1)),
         right_clip = as.integer(str_extract(cigar, ".*[^0-9]([0-9]+)S$", 1))) %>%
  replace_na(list(left_clip = 0, right_clip = 0)) %>%
  select(-cigar, -rnext, -pnext, -seq, -flag, -tlen, -qual, -score) %>%
  # remove duplicate genes
  ## anti_join(dups_df, by = "gene") %>%
  ## mutate(gene = if_else(is.na(new_gene), gene, new_gene)) %>%
  ## select(-new_gene) %>%
  # remove blacklisted genes
  ## filter(!gene %in% BLACKLIST) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  mutate(major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5)) %>%
  # get rid of stuff we probably don't care about
  filter((major == "V" & mapq > map_params$mapQ$V)
         | (major == "D" & mapq > map_params$mapQ$D)
         | (major == "J" & mapq > map_params$mapQ$J)
         | (major == "C" & mapq > map_params$mapQ$C)) %>%
  filter((major == "V" & scorefrac > map_params$score$V)
         | (major == "D" & scorefrac > map_params$score$D)
         | (major == "J" & scorefrac > map_params$score$J)
         | (major == "C" & minor != "DEL" & scorefrac > map_params$score$C)
         | (major == "C" & minor == "DEL" & scorefrac > 0.3)) %>%
  left_join(headers_df, by = c("gene", "allele")) %>%
  select(-gene) %>%
  # filter out genes that don't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  relocate(rname, pos, end) %>%
  # put all reads in the same orientation
  group_by(rname) %>%
  mutate(revfrac = mean(is_reversed),
         read_reversed = case_when(revfrac > 0.95 ~ TRUE,
                                   revfrac < 0.05 ~ FALSE,
                                   TRUE ~ NA),
         nV = sum(major == "V"),
         nD = sum(major == "D"),
         nJ = sum(major == "J"),
         nC = sum(major == "C"),
         .pos = if_else(read_reversed, rlen - end + 1, pos),
         .end = if_else(read_reversed, rlen - pos + 1, end),
         .right_clip = if_else(read_reversed, left_clip, right_clip),
         .left_clip = if_else(read_reversed, right_clip, left_clip),
         right_clip = .right_clip,
         left_clip = .left_clip,
         pos = .pos,
         end = .end) %>%
  select(-.end, -.pos, -.left_clip, -.right_clip) %>%
  # remove duplicated genes in the IGK locus
  filter(str_sub(minor, 2, 2) != "D")

ambiguous_orientation <- df %>%
  filter(is.na(pos))

merged_df <- df %>%
  filter(!is.na(pos)) %>%
  # do cheap bedtools merge thing; for every overlapping region, keep the
  # start/end of the first/last region
  arrange(rname, pos, end) %>%
  group_by(rname) %>%
  mutate(merge_pos = map2(pos, end, list) %>%
           accumulate(get_merge_end) %>%
           map_dbl(~ .x[[1]])) %>%
  group_by(rname, merge_pos) %>%
  mutate(merge_end = max(end)) %>%
  ungroup()

gene_pos_df <- merged_df %>%
  select(rname, merge_pos) %>%
  unique() %>%
  group_by(rname) %>%
  mutate(pos_index = row_number())

# TODO save this df
low_counts <- gene_pos_df %>%
  group_by(rname) %>%
  filter(n() < 3)

gene_df <- merged_df %>%
  anti_join(low_counts, by = "rname") %>%
  arrange(rname, merge_pos, gene_pos) %>%
  group_by(rname) %>%
  left_join(gene_pos_df, by = c("rname", "merge_pos"))

optimal_gene_df <- gene_df %>%
  # get rid of genes oriented against the majority
  filter(read_reversed == is_reversed) %>%
  group_by(rname, pos_index, gene_pos) %>%
  # for each gene at each location, only keep the highest scoring allele (this
  # will greatly reduce the number of combinations I need to deal with)
  filter(scorefrac == max(scorefrac)) %>%
  ungroup() %>%
  select(rname, pos_index, gene_pos, gene_index, scorefrac) %>%
  nest(reads = c(pos_index, gene_pos, gene_index, scorefrac)) %>%
  mutate(reads = map(reads,
                     optimize_order,
                     map_params$optimize$gap,
                     map_params$optimize$score)
         ) %>%
  unnest(reads) %>%
  left_join(gene_df, by = c("rname", "gene_index", "pos_index")) %>%
  relocate(rname, pos, end)

optimal_gene_df %>%
  group_by(rname) %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line() +
  labs(x = "read position", y = "gene rank", color = "region type")
ggsave(snakemake@output[["profile"]])

optimal_gene_df %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), group = rname, fill = major)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  labs(x = "gene rank", y = "read", fill = "region type")
ggsave(snakemake@output[["tiling"]])

vaf_df <- optimal_gene_df %>%
  mutate(gene = sprintf("%s%s", major, minor)) %>%
  group_by(gene, allele) %>%
  mutate(atot = n()) %>%
  group_by(gene) %>%
  mutate(vaf = atot / n()) %>%
  ungroup() %>%
  select(gene, allele, vaf) %>%
  arrange(gene, allele) %>%
  unique()

dedup_vaf_df <- optimal_gene_df %>%
  mutate(gene = sprintf("%s%s", major, minor)) %>%
  left_join(vaf_df, by = c("gene", "allele"))

dedup_vaf_df %>%
  select(rname, gene, gene_pos, allele, vaf) %>%
  group_by(gene) %>%
  mutate(allele = as.integer(fct_reorder(factor(allele), vaf, .desc = TRUE))) %>%
  ungroup() %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), fill = factor(allele))) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  labs(x = "gene rank", y = "read", fill = "allele index")
ggsave(snakemake@output[["tiling_allele"]])

dedup_vaf_df %>%
  write_tsv(snakemake@output[["tsv"]])

optimal_gene_df %>%
  mutate(.pos = if_else(read_reversed, rlen - end + 1, pos),
         .end = if_else(read_reversed, rlen - pos + 1, end),
         pos = .pos - 1,
         end = .end - 1,
         name = sprintf("%s%s*%s", major, minor, allele),
         score = scorefrac * 1000,
         strand = if_else(is_reversed, "-", "+")) %>%
  select(rname, pos, end, name, score, strand) %>%
  write_tsv(snakemake@output[["bed"]], col_names = FALSE)

ambiguous_orientation %>%
  select(-pos, -end, -left_clip, -right_clip, -read_reversed) %>%
  write_tsv(snakemake@output[["ambiguous"]])

low_counts %>%
  write_tsv(snakemake@output[["low_count"]])
  
