library(tidyverse)

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  rename(gene_pos = pos) %>%
  filter(gene != "IGHD") %>%
  filter(str_sub(gene, 4, 4) %in% c("V", "D", "J"))

headers_df <- read_tsv(
  "../../results/igmt/headers.tsv",
  col_names = c(
    "accession",
    "gene",
    "species",
    "functionality",
    "regiontype",
    "pos",
    "len",
    "codonstart",
    "nucaddstart",
    "nucaddend",
    "correction",
    "numAA",
    "numchars",
    "partial",
    "gene_revcomp"
  ),
  col_types = "c"
) %>%
  select(-accession, -species, -pos, -len, -codonstart, -correction,
         -numAA, -matches("nuc"), -partial, -numchars, -regiontype) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  full_join(order_df, by = "gene") %>%
  filter(!is.na(allele)) %>%
  filter(!is.na(gene_pos)) %>%
  filter(functionality != "ORF") %>%
  select(-functionality) %>%
  mutate(gene_revcomp = !is.na(gene_revcomp))

tags_df <- read_lines(
  "../../results/immuno_alignments/ont/split/IGH/all.tags"
) %>%
  tibble(tags = .) %>%
  mutate(id = row_number(), tags = str_split(tags, "\t")) %>%
  unnest(tags) %>%
  mutate(key = str_extract(tags, "^(.*):.*:.*", 1),
         value = str_extract(tags, ".*:.*:(.*)$", 1)) %>%
  pivot_wider(id_cols = id, names_from = key, values_from = value) %>%
  mutate(NM = as.numeric(NM),
         AS = as.numeric(AS),
         XS = as.numeric(XS)) %>%
  rename(edit_distance = NM,
         edit_string = MD,
         score = AS,
         subscore = XS) %>%
  select(-id)

sam_df <- readr::read_tsv(
  "../../results/immuno_alignments/ont/split/IGH/all.sam",
  col_names = c(
    "qname",
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
  mutate(gene = str_extract(qname, ".*\\|(.*)\\|.*", 1),
         revcomp = (flag %% 32) %/% 16 == 1,
         seqlen = str_length(seq),
         end = pos + seqlen,
         scorefrac = score / seqlen) %>%
  # get rid of stuff we probably don't care about
  filter(mapq > 30) %>%
  filter(scorefrac > 0.9) %>%
  select(-qname, -cigar, -rnext, -pnext, -edit_string, -seq, -flag, -tlen,
         -subscore, -qual, -score) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  mutate(major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5)) %>%
  left_join(headers_df, by = c("gene", "allele")) %>%
  select(-gene) %>%
  mutate(is_reversed = gene_revcomp != revcomp) %>%
  select(-matches("revcomp")) %>%
  # filter out genes that doesn't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  relocate(rname, pos, end) %>%
  # do cheap bedtools merge thing; for every overlapping region, keep the
  # start/end of the first/last region
  arrange(rname, pos, end) %>%
  group_by(rname) %>%
  mutate(.maxend = cummax(end),
         merge_pos = if_else(row_number() > 1 & pos < lag(.maxend), lag(pos), pos)) %>%
  group_by(rname, pos) %>%
  mutate(merge_end = max(end)) %>%
  select(-.maxend) %>%
  ungroup()

# find alleles that are mapped to multiple places; this should be almost or
# totally empty
dup_allele_df <- df %>%
  group_by(major, minor, allele, rname) %>%
  arrange(rname, major, minor, allele) %>%
  filter(n() > 1)

dedup_allele_df <- df %>%
  # if there are dup alleles, only keep ones with highest score
  group_by(major, minor, allele, rname) %>%
  filter(max(scorefrac) == scorefrac) %>%
  ungroup() %>%
  # for locations that have multiple genes, choose the one with highest score,
  # followed by min edit distance, followed by arbitrary first choice
  group_by(rname, merge_pos, major) %>%
  filter(max(scorefrac) == scorefrac) %>%
  filter(min(edit_distance) == edit_distance) %>%
  filter(row_number() == 1) %>%
  ungroup()

# find locations with multiple mappings; this should be empty
dup_location_df <- dedup_allele_df %>%
  group_by(rname, merge_pos) %>%
  filter(n() > 1)

## merge_df %>%
##   group_by(rname, merge_pos, major, minor) %>%
##   ## select(-functionality, -matches("merge"), -mapq)
##   filter(max(scorefrac) == scorefrac)

  ## relocate(rname, pos, end, mapq, scorefrac, edit_distance, seqlen) %>%

## perfect_df <- dedup_allele_df %>%
##   arrange(rname, pos, end) %>%
##   group_by(rname) %>%
##   mutate(rel_gene_pos = gene_pos - lag(gene_pos, default = 0)) %>%
##   filter(min(rel_gene_pos) > 0) %>%
##   ungroup()

## badV <- dedup_allele_df %>%
##   filter(major == "V") %>%
##   arrange(rname, pos, end) %>%
##   group_by(rname) %>%
##   mutate(rel_gene_pos = gene_pos - lag(gene_pos, default = 0)) %>%
##   filter(max(abs(rel_gene_pos)) > 50) %>%
##   ungroup() %>%
##   pull(rname) %>%
##   unique()

vdj_summary <- dedup_allele_df %>%
  group_by(rname) %>%
  summarize(nV = sum(major == "V"),
            nD = sum(major == "D"),
            nJ = sum(major == "J"))

hasJ_df <- dedup_allele_df %>%
  group_by(rname) %>%
  filter(n() > 2)
  ## filter(rname %in% (vdj_summary %>% filter(nJ > 1) %>% pull(rname)))

hasJ_df %>%
  mutate(revfrac = mean(is_reversed),
         revcat = case_when(revfrac > 0.95 ~ "reverse",
                            revfrac < 0.05 ~ "forward",
                            TRUE ~ "both")) %>%
  filter(revcat != "both") %>%
  ggplot(aes(gene_pos, merge_pos, group = rname, color = major)) +
  geom_point() +
  geom_line() +
  facet_wrap(c("revcat"))

hasJ_df %>%
  mutate(revfrac = mean(is_reversed),
         revcat = case_when(revfrac > 0.95 ~ "reverse",
                            revfrac < 0.05 ~ "forward",
                            TRUE ~ "both")) %>%
  filter(revcat != "both") %>%
  mutate(merge_pos = if_else(revcat == "forward", merge_pos, seqlen - merge_pos)) %>%
  ## group_by(rname) %>%
  ## arrange(merge_pos) %>%
  ## mutate(last_gene_pos = lag(gene_pos),
  ##        gap = if_else(is.na(last_gene_pos), 0, gene_pos - last_gene_pos - 1)) %>%
  ## ungroup() %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), group = rname, fill = major)) +
  geom_tile() +
  theme(axis.text.y = element_blank())

## hasJ_df %>%
##   group_by(gene_pos, major) %>%
##   tally() %>%
##   ggplot(aes(gene_pos, n, fill = major)) +
##   geom_col()

# make bed file
hasJ_df %>%
  mutate(gene = sprintf("%s%s*%s", major, minor, allele)) %>%
  select(-major, -minor, -allele, -pos, -end, -seqlen) %>%
  relocate(rname, merge_pos, merge_end, gene)
