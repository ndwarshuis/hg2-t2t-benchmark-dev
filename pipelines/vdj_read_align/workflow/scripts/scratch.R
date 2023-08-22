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
         | (major != "D" & ((mapq > 30 & !is_hifi) | (mapq > 40 & is_hifi)))) %>%
  filter((major == "V" & ((scorefrac > 0.85 & !is_hifi) | (scorefrac > 0.95 & is_hifi)))
         | (major == "D" & ((scorefrac > 0.55) & !is_hifi) | (scorefrac > 0.85 & is_hifi))
         | (major == "J" & ((scorefrac > 0.5) & !is_hifi) | (scorefrac > 0.80 & is_hifi))
         | (major == "C" & minor != "DEL" & scorefrac > 0.5)
         | (major == "C" & minor == "DEL" & scorefrac > 0.3)
         ) %>%
  ## left_join(headers_df, by = c("gene", "allele")) %>%
  left_join(order_df, by = c("gene")) %>%
  select(-gene) %>%
  ## mutate(is_reversed = gene_revcomp != revcomp) %>%
  ## select(-matches("revcomp")) %>%
  # filter out genes that don't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  relocate(rname, pos, end) %>%
  # put all reads in the same orientation
  group_by(rname) %>%
  mutate(revfrac = mean(is_reversed),
         revcat = case_when(revfrac > 0.95 ~ "reverse",
                            revfrac < 0.05 ~ "forward",
                            TRUE ~ "both")) %>%
  filter(revcat != "both") %>%
  # get rid of genes oriented against the majority
  filter((revcat == "forward" & !is_reversed)
         | (revcat == "reverse" & is_reversed)) %>%
  mutate(.pos = if_else(revcat == "reverse", rlen - end + 1, pos),
         .end = if_else(revcat == "reverse", rlen - pos + 1, end),
         .right_clip = if_else(revcat == "reverse", left_clip, right_clip),
         .left_clip = if_else(revcat == "reverse", right_clip, left_clip),
         right_clip = .right_clip,
         left_clip = .left_clip,
         pos = .pos,
         end = .end) %>%
  select(-.end, -.pos, -.left_clip, -.right_clip) %>%
  # do cheap bedtools merge thing; for every overlapping region, keep the
  # start/end of the first/last region
  arrange(rname, pos, end) %>%
  group_by(rname) %>%
  mutate(merge_pos = map2(pos, end, list) %>%
           accumulate(get_merge_end) %>%
           map_dbl(~ .x[[1]])) %>%
  group_by(rname, merge_pos) %>%
  mutate(merge_end = max(end)) %>%
  ungroup() %>%
  # remove duplicated genes in the IGK locus
  filter(str_sub(minor, 2, 2) != "D")

gene_df <- df %>%
  arrange(rname, merge_pos, gene_pos) %>%
  group_by(rname) %>%
  mutate(gene = paste0(major, minor)) %>%
  select(rname, merge_pos, gene_pos) %>%
  unique()

optimize_order <- function(df) {
  xs <- df %>%
    group_by(pos_index) %>%
    group_map(~ .x$gene_pos)
  branches <- accumulate(xs, ~ .x * length(.y), .init = 1)[-1]
  ngenes <- length(xs)
  ncombos <- last(branches)
  paths <- xs %>%
    map2(branches, ~ rep(.x, .y / length(.x), each = ncombos / .y)) %>%
    flatten_dbl() %>%
    matrix(nrow = ncombos)
  dups <- apply(paths, 1, function(x) any(duplicated(x)))
  unique_paths <- paths[!dups, ]
  best_path <- if (is.null(dim(unique_paths))) {
    unique_paths
  } else {
    ncombos_nodup <- nrow(unique_paths)
    gaps <- (unique_paths[, -ncol(unique_paths)] - unique_paths[, -1]) %>%
      abs() %>%
      magrittr::add(-1) %>%
      matrix(nrow = ncombos_nodup)
    scores <- apply(gaps, 1, sum)
    best_scores <- scores == min(scores)
    best_paths <- unique_paths[best_scores, ] %>%
      matrix(nrow = sum(best_scores))
    best_paths[1, ]
  }
  tibble(gene_pos = best_path) %>%
    mutate(pos_index = row_number())
}

optimal_gene_df <- gene_df %>%
  select(rname, merge_pos) %>%
  unique() %>%
  group_by(rname) %>%
  mutate(pos_index = row_number()) %>%
  # remove reads with only 1 gene on them
  filter(n() > 1) %>%
  ## filter(rname == "@0a1a25e5-0f65-4f42-9a0a-cb48382e3c84") %>%
  left_join(gene_df, by = c("rname", "merge_pos")) %>%
  select(-merge_pos) %>%
  nest(reads = c(pos_index, gene_pos)) %>%
  mutate(reads = map(reads, optimize_order)) %>%
  unnest(reads)
  ## pull(reads) %>%
  ## first() %>%
  ## group_by(pos_index) %>%
  ## group_map(~ .x$gene_pos)

gene_df %>%
  mutate(pos_index = as.integer(factor(merge_pos))) %>%
  left_join(optimal_gene_df, by = c("rname", "gene_pos")) %>%
  View()

optimal_gene_df %>%
  select(pos_index, gene_pos)

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

dist_df <- df %>%
  pivot_longer(cols = -gene, values_to = "distance", names_to = "other") %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = c(gene, other), values_to = "name", names_to = "var") %>%
  mutate(var = if_else(var == "gene", "geneA", "geneB")) %>%
  mutate(major = str_sub(name, 4, 4),
         minor = str_sub(str_extract(name, "(.*)\\*.*", 1), 5),
         allele = as.integer(str_extract(name, ".*\\*(.*)", 1))) %>%
  filter(str_extract(name, "(.*)\\*.*", 1) %in% genes_with_rank) %>%
  pivot_wider(id_cols = c(id, distance),
              names_from = var,
              values_from = c(name, major, minor, allele),
              names_glue = "{var}_{.value}") %>%
  filter(!is.na(geneA_name)) %>%
  filter(!is.na(geneB_name)) %>%
  select(-id)

dist_df %>%
  filter(geneA_minor != geneB_minor) %>%
  filter(distance < 10) %>%
  ggplot(aes(geneA_name, geneB_name, fill = distance)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

sim_df <- dist_df %>%
  filter(geneA_minor != geneB_minor) %>%
  filter(distance < 10) %>%
  arrange(geneA_name) %>%
  mutate(rankA = as.numeric(factor(geneA_name))) %>%
  arrange(geneB_name) %>%
  mutate(rankB = as.numeric(factor(geneB_name))) %>%
  filter(rankA < rankB) %>%
  select(geneA_name, geneB_name, rankA, rankB)

similar_genes <- sim_df %>%
  anti_join(select(sim_df, rankB) %>% unique(), by = c("rankA" = "rankB")) %>%
  group_by(geneA_name) %>%
  group_map(~ c(.y$geneA_name, .x$geneB_name))

partition <- function(xs, fun) {
  # gross...
  list(pass = keep(xs, fun), fail = discard(xs, fun))
}

pick_next <- function(xs, fun, more_funs = list()) {
  pick_first <- function(xs) {
    ys <- sort(xs)
    list(rep = ys[[1]], drop = ys[-1])
  }
  next_fun <- if (length(more_funs) == 0) {
    pick_first
  } else {
    function(ys) pick_next(ys, more_funs[[1]], more_funs[-1])
  }
  res <- partition(xs, fun)
  if (length(res$pass) == 1) {
    list(rep = res$pass[[1]], drop = res$fail)
  } else if (length(res$pass) == 0) {
    next_fun(xs)
  } else {
    .res <- next_fun(res$pass)
    .res$drop = c(.res$drop, res$fail)
    .res
  }
}

is_not_dup <- function(gene) {
  str_sub(str_extract(gene, "(.*)\\*.*", 1), -1) != "D"
}

is_not_K_dup <- function(gene) {
  !str_detect(str_sub(gene, 4), "D-")
}

similar_genes %>%
  map(~ pick_next(.x, ~ !str_detect(.x, "OR"), list(is_not_dup, is_not_K_dup))) %>%
  imap_dfr(~ tibble(id = .y, new = .x$rep, old = .x$drop)) %>%
  View()
