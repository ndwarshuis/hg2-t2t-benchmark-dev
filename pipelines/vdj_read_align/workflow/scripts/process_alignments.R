library(tidyverse)

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  rename(gene_pos = pos)

headers_df <- read_tsv(
  "../../results/igmt/IGH_headers.fa",
  col_names = c("gene", "allele", "functionality", "gene_orientation"),
  col_types = "c"
) %>%
  inner_join(order_df, by = "gene") %>%
  filter(functionality != "ORF") %>%
  select(-functionality) %>%
  mutate(gene_revcomp = gene_orientation == "-") %>%
  select(-gene_orientation)

tags_df <- read_lines(
  "../../results/immuno_alignments/ont/split/IGH/all_allele.tags"
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

df <- sam_df %>%
  mutate(revcomp = (flag %% 32) %/% 16 == 1,
         seqlen = str_length(seq),
         end = pos + seqlen,
         scorefrac = score / seqlen) %>%
  # get rid of stuff we probably don't care about
  filter(mapq > 30) %>%
  filter(scorefrac > 0.9) %>%
  select(-cigar, -rnext, -pnext, -edit_string, -seq, -flag, -tlen,
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
            nJ = sum(major == "J"),
            nC = sum(major == "C"))

VJC_df <- vdj_summary %>%
  filter(nV > 1, nJ > 1, nC > 1)

hasJ_df <- dedup_allele_df %>%
  group_by(rname) %>%
  filter(n() > 2)
  ## filter(rname %in% (vdj_summary %>% filter(nJ > 1) %>% pull(rname)))

hasJ_df %>%
  filter(rname %in% (VJC_df %>% pull(rname))) %>%
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

# figure out which alleles are likely wrong

ngene_df <- dedup_allele_df %>%
  group_by(major, minor, gene_pos) %>%
  summarize(ngene = n(), .groups = "drop") %>%
  arrange(desc(ngene))

nallele_df <- dedup_allele_df %>%
  group_by(major, minor, allele) %>%
  summarize(nallele = n(), .groups = "drop") %>%
  left_join(ngene_df, by = c("major", "minor")) %>%
  mutate(vaf = nallele / ngene) %>%
  arrange(major, minor, vaf)

nallele_df %>%
  filter(vaf < 1) %>%
  mutate(gene = sprintf("%s%s", major, minor)) %>%
  ggplot(aes(vaf, gene, fill = allele)) +
  geom_col()

flank_df <- readr::read_tsv(
  "../../results/immuno_alignments/ont/split/IGH/all_motif.sam",
  col_types = "cicidcciicccccc",
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
    "qual",
    "XA",
    "MD",
    "NM",
    "XM"
  )) %>%
  mutate(edit = as.numeric(str_replace(NM, "NM:i:", ""))) %>%
  select(qname, flag, rname, pos, edit)

get_rss <- function(df, forward) {
  # filter by forward/reverse
  flag_df <- df %>%
    filter(flag == ifelse(forward, 0, 16)) %>%
    select(-flag)
  # grab all the 5p motifs, and make columns of positions within which the other
  # motif should fall if it is valid
  offset <- ifelse(forward, 7, 9)
  df5p <- flag_df %>%
    filter(qname == ifelse(forward, "heptamer", "nonamer")) %>%
    mutate(lower12 = pos + offset + 10,
           upper12 = pos + offset + 14,
           lower23 = pos + offset + 21,
           upper23 = pos + offset + 25) %>%
    select(-qname)
  # grab all 3p motifs and relabel columns so they don't clash
  df3p <- flag_df %>%
    filter(qname == ifelse(forward, "nonamer", "heptamer")) %>%
    select(-qname) %>%
    rename(pos3p = pos, edit3p = edit)
  # do a fuzzy inner join on the 12bp positions
  df12 <- df5p %>%
    inner_join(df3p, join_by(rname, between(y$pos3p, x$lower12, x$upper12))) %>%
    mutate(is12 = TRUE)
  # ...and same with the 23bp positions
  df23 <- df5p %>%
    inner_join(df3p, join_by(rname, between(y$pos3p, x$lower23, x$upper23))) %>%
    mutate(is12 = FALSE)
  # combine rows and summarize
  bind_rows(df12, df23) %>%
    mutate(
      edit = edit + edit3p,
      strand = if_else(forward, "+", "-"),
      spacer = pos3p - pos - offset,
      end = pos3p + ifelse(forward, 9, 7),
      thickStart = pos + ifelse(forward, 7, 9),
      thickEnd = end - ifelse(forward, 9, 7)
    ) %>%
  select(-matches("lower"), -matches("upper"), -matches("3p"))
}

rss_df <- get_rss(flank_df, TRUE) %>%
  bind_rows(get_rss(flank_df, FALSE)) %>%
  arrange(rname, pos) %>%
  rename("#rname" = rname,
         start = pos) %>%
  mutate(name = if_else(is12, "RSS12", "RSS23"),
         score = 1000 - edit * 200) %>%
  relocate(`#rname`, start, end, name, score, strand, thickStart, thickEnd)
         
