suppressMessages(library(tidyverse))

order_df <- read_tsv(
  snakemake@input[["order"]],
  col_types = "ci"
) %>%
  rename(gene_pos = pos,
         name = gene)

headers_df <- read_tsv(
  snakemake@input[["headers"]],
  col_names = c("name", "allele", "functionality", "gene_orientation"),
  col_types = "c"
) %>%
  inner_join(order_df, by = "name") %>%
  filter(functionality != "ORF") %>%
  select(-functionality) %>%
  mutate(gene_revcomp = gene_orientation == "-") %>%
  select(-gene_orientation)

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
  rename(edit = NM,
         score = AS) %>%
  select(-id)

sam_df <- readr::read_tsv(
  snakemake@input[["sam"]],
  col_names = c(
    "qname",
    "flag",
    "#rname",
    "start",
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

sam_df %>%
  mutate(name = str_extract(qname, ".*\\|(.*)\\|.*", 1),
         revcomp = (flag %% 32) %/% 16 == 1,
         strand = if_else(revcomp, "-", "+"),
         end = start + str_length(seq),
         scorefrac = score / str_length(seq),
         score = 1000 * scorefrac,
         thickStart = start,
         thickEnd = end) %>%
  # get rid of stuff we probably don't care about
  filter(mapq > 30) %>%
  filter(scorefrac > 0.9) %>%
  select(-qname, -cigar, -rnext, -pnext, -seq, -flag, -tlen, -qual) %>%
  left_join(headers_df, by = "name") %>%
  mutate(is_reversed = gene_revcomp != revcomp) %>%
  select(-revcomp) %>%
  # filter out genes that doesn't have a known order (not as useful)
  filter(!is.na(gene_pos)) %>%
  # do cheap bedtools merge thing; for every overlapping region, keep the
  # start/end of the first/last region
  ## arrange(rname, pos, end) %>%
  ## group_by(rname) %>%
  ## mutate(.maxend = cummax(end),
  ##        merge_pos = if_else(row_number() > 1 & pos < lag(.maxend), lag(pos), pos)) %>%
  ## group_by(rname, pos) %>%
  ## mutate(merge_end = max(end)) %>%
  ## select(-.maxend) %>%
  ## ungroup() %>%
  relocate(`#rname`, start, end, name, score, thickStart, thickEnd) %>%
  write_tsv(snakemake@output[[1]])
