library(tidyverse)

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  rename(gene_pos = pos)

constants <- c(
  paste0("IGHC", c(paste0("G", c(1:4, "P")), "E", "M", "D", "A1", "A2", "EP1", "EP2")),
  "IGKC",
  paste0("IGLC", 1:7),
  "TRAC",
  "TRDC",
  "TRBC1",
  "TRBC2",
  "TRGC1",
  "TRGC2"
) %>% tibble(gene = ., allele = "01", functionality = "", gene_orientation = "+")

headers_df <- read_tsv(
  "../../results/igmt/IGK_vdj_headers.fa",
  col_names = c("gene", "allele", "functionality", "gene_orientation"),
  col_types = "c"
) %>%
  bind_rows(constants) %>%
  inner_join(order_df, by = "gene") %>%
  filter(functionality != "ORF") %>%
  select(-functionality) %>%
  mutate(gene_revcomp = gene_orientation == "-") %>%
  select(-gene_orientation)

reads_df <- read_tsv(
  "../../results/immuno_alignments/ont/split/IGK.fastq.fai",
  col_names = FALSE,
  col_types = c(X1 = "c", X2 = "i", .default = "-")) %>%
  rename(rname = X1, rlen = X2) %>%
  mutate(rname = paste0("@", rname))

tags_df <- read_lines(
  "../../results/immuno_alignments/ont/split/IGK/all_allele.tags"
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
  "../../results/immuno_alignments/ont/split/IGK/all_allele.sam",
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

get_merge_end <- function(last, x) {
  if (x[[1]] > last[[2]]) {
    x
  } else {
    last
  }
}

df <- sam_df %>%
  left_join(reads_df, by = "rname") %>%
  mutate(revcomp = (flag %% 32) %/% 16 == 1,
         seqlen = str_length(seq),
         end = map2_dbl(pos + seqlen, rlen, min),
         scorefrac = score / seqlen) %>%
  select(-cigar, -rnext, -pnext, -edit_string, -seq, -flag, -tlen,
         -subscore, -qual, -score) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  mutate(major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5),
         rss = if_else(major == "C", TRUE, str_detect(gene, "-RSS")),
         gene = if_else(major == "C", gene, str_extract(gene, "(.*)-(-no)?RSS", 1))) %>%
  # get rid of stuff we probably don't care about
  filter(mapq > 30) %>%
  filter((major == "V" & scorefrac > 0.9)
         | (major == "D" & scorefrac > 0.65)
         | (major == "J" & scorefrac > 0.75)
         | (major == "C" & scorefrac > 0.5)
         ) %>%
  left_join(headers_df, by = c("gene", "allele")) %>%
  select(-gene) %>%
  mutate(is_reversed = gene_revcomp != revcomp) %>%
  select(-matches("revcomp")) %>%
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
  mutate(.pos = if_else(revcat == "reverse", rlen - end + 1, pos),
         .end = if_else(revcat == "reverse", rlen - pos + 1, end),
         pos = .pos,
         end = .end) %>%
  select(-.end, -.pos) %>%
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
  filter(str_sub(minor, 2, 2) != "D")

# find alleles that are mapped to multiple places; this should be almost or
# totally empty
dup_allele_df <- df %>%
  group_by(major, minor, allele, rname) %>%
  arrange(rname, major, minor, allele) %>%
  filter(n() > 1)

dedup_allele_df <- df %>%
  filter(rss) %>%
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
  ungroup() %>%
  group_by(rname) %>%
  mutate(nV = sum(major == "V"),
         nD = sum(major == "D"),
         nJ = sum(major == "J"),
         nC = sum(major == "C")) %>%
  ungroup()

# find locations with multiple mappings; this should be empty
dup_location_df <- dedup_allele_df %>%
  group_by(rname, merge_pos) %>%
  filter(n() > 1)

lastV <- dedup_allele_df %>%
  group_by(rname) %>%
  filter(major == "V") %>%
  filter(max(pos) == pos) %>%
  mutate(Vname = paste0(major, minor)) %>%
  select(rname, Vname)

firstJ <- dedup_allele_df %>%
  group_by(rname) %>%
  filter(major == "J") %>%
  filter(min(pos) == pos) %>%
  mutate(Jname = paste0(major, minor)) %>%
  select(rname, Jname)

firstC <- dedup_allele_df %>%
  group_by(rname) %>%
  filter(major == "C") %>%
  filter(min(pos) == pos) %>%
  mutate(Cname = paste0(major, minor)) %>%
  select(rname, Cname)

VJC_combos <- lastV %>%
  inner_join(firstJ, by = "rname") %>%
  inner_join(firstC, by = "rname") %>%
  arrange(Vname, Jname, Cname)

VJC_combos %>%
  mutate(combo = str_c(Vname, Jname, Cname, sep = "_")) %>%
  ggplot(aes(y = combo)) +
  geom_bar()

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

dedup_allele_df %>%
  group_by(rname) %>%
  ## filter(min(gene_pos) < 22 & max(gene_pos) > 44) %>%
  ## filter(gene_pos > 60) %>%
  ## filter(rname %in% (vdj_summary %>% filter(nJ > 1) %>% pull(rname))) %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()
  ## facet_wrap(c("rname"))

dedup_allele_df %>%
  filter(rname == "@187b9e81-d9a1-485a-9e1a-a7157df7ea5c")

## dedup_allele_df %>%
##   group_by(rname) %>%
##   filter(min(gene_pos) < 22 & max(gene_pos) > 44) %>%
##   ungroup() %>%
##   arrange(rname, gene_pos) %>%
##   select(rname)
##   ## ggplot(aes(gene_pos, pos, group = rname, color = major)) +
##   ## geom_point() +
##   ## geom_line()
  

hasJ_df %>%
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
