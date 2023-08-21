suppressMessages(library(tidyverse))

BLACKLIST <- c("IGHV(II)-30-21-RSS*01")

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

get_merge_end <- function(last, x) {
  if (x[[1]] > last[[2]]) {
    x
  } else {
    last
  }
}

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
         seqlen = str_length(seq),
         end = map2_dbl(pos + seqlen, rlen, min),
         scorefrac = score / seqlen) %>%
  select(-cigar, -rnext, -pnext, -seq, -flag, -tlen, -qual, -score) %>%
  # remove duplicate genes
  left_join(dups_df, by = "gene") %>%
  mutate(gene = if_else(is.na(new_gene), gene, new_gene)) %>%
  select(-new_gene) %>%
  # remove blacklisted genes
  ## filter(!gene %in% BLACKLIST) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  mutate(major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5),
         rss = if_else(major == "C", TRUE, str_detect(gene, "-RSS")),
         gene = if_else(major == "C", gene, str_extract(gene, "(.*)-(-no)?RSS", 1))) %>%
  # get rid of stuff we probably don't care about
  filter((major == "D" & mapq > 10)
         | (major != "D" & mapq > 30)) %>%
  filter((major == "V" & scorefrac > 0.85)
         | (major == "D" & scorefrac > 0.55)
         | (major == "J" & scorefrac > 0.5)
         | (major == "C" & minor != "DEL" & scorefrac > 0.5)
         | (major == "C" & minor == "DEL" & scorefrac > 0.3)
         ) %>%
  left_join(headers_df, by = c("gene", "allele")) %>%
  select(-gene) %>%
  ## mutate(is_reversed = gene_revcomp != revcomp) %>%
  ## select(-matches("revcomp")) %>%
  # filter out genes that don't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  relocate(rname, pos, end) %>%
  # put all reads in the same orientation
  group_by(rname) %>%
  mutate(revfrac = mean(is_reversed),
         revcat = case_when(revfrac > 0.90 ~ "reverse",
                            revfrac < 0.10 ~ "forward",
                            TRUE ~ "both")) %>%
  filter(revcat != "both") %>%
  filter((revcat == "forward" & !is_reversed)
         | (revcat == "reverse" & is_reversed)) %>%
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
  # remove duplicated genes in the IGK locus
  filter(str_sub(minor, 2, 2) != "D")

dedup_df <- df %>%
  filter(rss) %>%
  # if there are dup alleles, only keep ones with highest score
  group_by(major, minor, rname) %>%
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
  filter(n() > 2) %>%
  ungroup() %>%
  select(-rss) %>%
  mutate(minor = str_replace(minor, "-RSS", ""))

dedup_df %>%
  group_by(rname) %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line() +
  labs(x = "read position", y = "gene rank", color = "region type")
ggsave(snakemake@output[["profile"]])

dedup_df %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), group = rname, fill = major)) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  labs(x = "gene rank", y = "read", fill = "region type")
ggsave(snakemake@output[["tiling"]])

vaf_df <- dedup_df %>%
  filter(major != "C") %>%
  mutate(gene = sprintf("%s%s", major, minor)) %>%
  group_by(gene, allele) %>%
  mutate(atot = n()) %>%
  group_by(gene) %>%
  mutate(vaf = atot / n()) %>%
  ungroup() %>%
  select(gene, allele, vaf) %>%
  arrange(gene, allele) %>%
  unique()

dedup_df %>%
  filter(major != "C") %>%
  mutate(gene = sprintf("%s%s", major, minor)) %>%
  select(rname, gene, gene_pos, allele) %>%
  left_join(vaf_df, by = c("gene", "allele")) %>%
  group_by(gene) %>%
  mutate(allele = as.integer(fct_reorder(factor(allele), vaf, .desc = TRUE))) %>%
  ungroup() %>%
  ggplot(aes(gene_pos, fct_reorder(rname, gene_pos), fill = factor(allele))) +
  geom_tile() +
  theme(axis.text.y = element_blank()) +
  labs(x = "gene rank", y = "read", fill = "allele index")
ggsave(snakemake@output[["tiling_allele"]])

dedup_df %>%
  write_tsv(snakemake@output[["tsv"]])

dedup_df %>%
  mutate(.pos = if_else(revcat == "reverse", rlen - end + 1, pos),
         .end = if_else(revcat == "reverse", rlen - pos + 1, end),
         pos = .pos,
         end = .end,
         name = sprintf("%s%s*%s", major, minor, allele),
         score = scorefrac * 1000,
         strand = if_else(is_reversed, "-", "+")) %>%
  select(rname, pos, end, name, score, strand) %>%
  write_tsv(snakemake@output[["bed"]], col_names = FALSE)
