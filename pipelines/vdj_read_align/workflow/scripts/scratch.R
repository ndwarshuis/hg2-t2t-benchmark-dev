library(tidyverse)

order_df <- read_tsv("../../static/gene_order.tsv") %>%
  rename(gene_pos = pos)

major_bounds <- order_df %>%
  mutate(major = str_sub(gene, 4, 4),
         locus = str_sub(gene, 1, 3)) %>%
  filter(major == "V" | major == "J") %>%
  group_by(major, locus) %>%
  mutate(maxV = max(gene_pos))

df <- readr::read_tsv(
  "../../results/immuno_alignments/ont/split/IGH/final/allele.tsv"
)

df %>%
  filter(major != "C") %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), group = rname, fill = left_clip)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  labs(x = "gene rank", y = "read", fill = "region type")

df %>%
  filter(major == "V") %>%
  filter(left_clip > 0) %>%
  select(rname, pos, gene, left_clip)

reads_df <- read_tsv(
  "../../results/immuno_alignments/ont/split/IGH.fastq.fai",
  col_names = FALSE,
  col_types = c(X1 = "c", X2 = "i", .default = "-")
) %>%
  rename(rname = X1, rlen = X2) %>%
  mutate(rname = paste0("@", rname))

tags_df <- read_lines(
  "../../results/immuno_alignments/ont/split/IGH/all_allele.tags",
) %>%
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
  "../../results/immuno_alignments/ont/split/IGH/all_allele.sam",
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

is_hifi <- FALSE

get_merge_end <- function(last, x) {
  if (x[[1]] > last[[2]]) {
    x
  } else {
    last
  }
}

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
  filter((major == "D" & ((mapq > 10 & !is_hifi) | (mapq > 30 & is_hifi)))
         | (major != "D" & ((mapq > 15 & !is_hifi) | (mapq > 40 & is_hifi)))) %>%
  filter((major == "V" & ((scorefrac > 0.85 & !is_hifi) | (scorefrac > 0.95 & is_hifi)))
         | (major == "D" & ((scorefrac > 0.55) & !is_hifi) | (scorefrac > 0.85 & is_hifi))
         | (major == "J" & ((scorefrac > 0.5) & !is_hifi) | (scorefrac > 0.80 & is_hifi))
         | (major == "C" & minor != "DEL" & scorefrac > 0.5)
         | (major == "C" & minor == "DEL" & scorefrac > 0.3)
         ) %>%
  ## left_join(headers_df, by = c("gene", "allele")) %>%
  left_join(order_df, by = c("gene")) %>%
  select(-gene) %>%
  # filter out genes that don't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  relocate(rname, pos, end) %>%
  # put all reads in the same orientation
  group_by(rname) %>%
  mutate(revfrac = mean(is_reversed),
         read_reversed = case_when(revfrac > 0.95 ~ TRUE,
                                   revfrac < 0.05 ~ FALSE,
                                   TRUE ~ NA)) %>%
  mutate(.pos = if_else(read_reversed, rlen - end + 1, pos),
         .end = if_else(read_reversed, rlen - pos + 1, end),
         .right_clip = if_else(read_reversed, left_clip, right_clip),
         .left_clip = if_else(read_reversed, right_clip, left_clip),
         right_clip = .right_clip,
         left_clip = .left_clip,
         pos = .pos,
         end = .end) %>%
  select(-.end, -.pos, -.left_clip, -.right_clip)
  # remove duplicated genes in the IGK locus
  ## filter(str_sub(minor, 2, 2) != "D")

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

low_counts <- gene_pos_df %>%
  group_by(rname) %>%
  filter(n() < 3)

gene_df <- merged_df %>%
  anti_join(low_counts, by = "rname") %>%
  arrange(rname, merge_pos, gene_pos) %>%
  group_by(rname) %>%
  left_join(gene_pos_df, by = c("rname", "merge_pos"))
  ## select(rname, pos_index, major, minor, allele, gene_index, gene_pos, scorefrac) %>%
  ## unique()

optimize_order <- function(df) {
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
  scores <- gap_scores * 0.5 + align_scores * 0.5
  best_combo <- detect_index(scores == min(scores), ~ .x == TRUE)
  tibble(gene_index = genes[best_combo, ]) %>%
    mutate(pos_index = as.integer(names(xs)))
}

optimal_gene_df <- gene_df %>%
  # get rid of genes oriented against the majority
  filter(read_reversed == is_reversed) %>%
  group_by(rname, pos_index, gene_pos) %>%
  # for each gene at each location, only keep the highest scoring allele (this
  # will greatly reduce the number of combinations I need to deal with)
  filter(scorefrac == max(scorefrac)) %>%
  ungroup() %>%
  ## filter(rname == "@665b90f0-8597-431d-b968-5687965c7be8") %>%
  select(rname, pos_index, gene_pos, gene_index, scorefrac) %>%
  nest(reads = c(pos_index, gene_pos, gene_index, scorefrac)) %>%
  mutate(reads = map(reads, optimize_order)) %>%
  unnest(reads) %>%
  left_join(gene_df, by = c("rname", "gene_index", "pos_index")) %>%
  relocate(rname, pos, end)

optimal_gene_df %>%
  ggplot(aes(pos, gene_pos, color = major, group = rname)) +
  geom_point() +
  geom_line()

## vaf_df <- df %>%
##   mutate(gene = sprintf("%s%s", major, minor)) %>%
##   group_by(gene, allele) %>%
##   mutate(atot = n()) %>%
##   group_by(gene) %>%
##   mutate(vaf = atot / n()) %>%
##   ungroup() %>%
##   select(gene, allele, vaf) %>%
##   arrange(gene, allele) %>%
##   unique()

## has110 <- df %>%
##   filter(gene_pos == 35) %>%
##   mutate(allele110 = as.numeric(factor(allele))) %>%
##   select(rname, allele110)

## df %>%
##   filter(major != "C") %>%
##   right_join(has110, by = "rname") %>%
##   group_by(gene) %>%
##   mutate(allele = as.integer(fct_reorder(factor(allele), vaf, .desc = TRUE))) %>%
##   ungroup() %>%
##   ggplot(aes(gene_pos, fct_reorder(rname, allele110), fill = factor(allele))) +
##   geom_tile() +
##   theme(axis.text.y = element_blank()) +
##   labs(x = "gene rank", y = "read", fill = "allele index")

## df %>%
##   filter(major != "C") %>%
##   mutate(gene = sprintf("%s%s", major, minor)) %>%
##   group_by(gene, allele) %>%
##   mutate(atot = n()) %>%
##   group_by(gene) %>%
##   mutate(vaf = atot / n()) %>%
##   ungroup() %>%
##   select(gene, allele, vaf) %>%
##   arrange(gene, allele) %>%
##   unique() %>%
##   group_by(gene) %>%
##   mutate(nallele = n()) %>%
##   filter(nallele > 1) %>%
##   ggplot(aes(vaf, gene, fill = allele)) +
##   geom_col(position = "stack") +
##   facet_wrap(c("nallele"), scale = "free_y")

## dist_df <- df %>%
##   pivot_longer(cols = -gene, values_to = "distance", names_to = "other") %>%
##   mutate(id = row_number()) %>%
##   pivot_longer(cols = c(gene, other), values_to = "name", names_to = "var") %>%
##   mutate(var = if_else(var == "gene", "geneA", "geneB")) %>%
##   mutate(major = str_sub(name, 4, 4),
##          minor = str_sub(str_extract(name, "(.*)\\*.*", 1), 5),
##          allele = as.integer(str_extract(name, ".*\\*(.*)", 1))) %>%
##   filter(str_extract(name, "(.*)\\*.*", 1) %in% genes_with_rank) %>%
##   pivot_wider(id_cols = c(id, distance),
##               names_from = var,
##               values_from = c(name, major, minor, allele),
##               names_glue = "{var}_{.value}") %>%
##   filter(!is.na(geneA_name)) %>%
##   filter(!is.na(geneB_name)) %>%
##   select(-id)

## dist_df %>%
##   filter(geneA_minor != geneB_minor) %>%
##   filter(distance < 10) %>%
##   ggplot(aes(geneA_name, geneB_name, fill = distance)) +
##   geom_tile() +
##   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## sim_df <- dist_df %>%
##   filter(geneA_minor != geneB_minor) %>%
##   filter(distance < 10) %>%
##   arrange(geneA_name) %>%
##   mutate(rankA = as.numeric(factor(geneA_name))) %>%
##   arrange(geneB_name) %>%
##   mutate(rankB = as.numeric(factor(geneB_name))) %>%
##   filter(rankA < rankB) %>%
##   select(geneA_name, geneB_name, rankA, rankB)

## similar_genes <- sim_df %>%
##   anti_join(select(sim_df, rankB) %>% unique(), by = c("rankA" = "rankB")) %>%
##   group_by(geneA_name) %>%
##   group_map(~ c(.y$geneA_name, .x$geneB_name))

## partition <- function(xs, fun) {
##   # gross...
##   list(pass = keep(xs, fun), fail = discard(xs, fun))
## }

## pick_next <- function(xs, fun, more_funs = list()) {
##   pick_first <- function(xs) {
##     ys <- sort(xs)
##     list(rep = ys[[1]], drop = ys[-1])
##   }
##   next_fun <- if (length(more_funs) == 0) {
##     pick_first
##   } else {
##     function(ys) pick_next(ys, more_funs[[1]], more_funs[-1])
##   }
##   res <- partition(xs, fun)
##   if (length(res$pass) == 1) {
##     list(rep = res$pass[[1]], drop = res$fail)
##   } else if (length(res$pass) == 0) {
##     next_fun(xs)
##   } else {
##     .res <- next_fun(res$pass)
##     .res$drop = c(.res$drop, res$fail)
##     .res
##   }
## }

## is_not_dup <- function(gene) {
##   str_sub(str_extract(gene, "(.*)\\*.*", 1), -1) != "D"
## }

## is_not_K_dup <- function(gene) {
##   !str_detect(str_sub(gene, 4), "D-")
## }

## similar_genes %>%
##   map(~ pick_next(.x, ~ !str_detect(.x, "OR"), list(is_not_dup, is_not_K_dup))) %>%
##   imap_dfr(~ tibble(id = .y, new = .x$rep, old = .x$drop)) %>%
##   View()
