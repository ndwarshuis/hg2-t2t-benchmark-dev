suppressMessages(library(tidyverse))

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  rename(gene_pos = pos)

df <- read_lines("../../resources/imgt_human.fa") %>%
  keep(~ str_starts(.x, ">")) %>%
  map(~ str_replace(.x, "^>", "")) %>%
  tibble(X1 = .) %>%
  separate(X1, sep = "\\|", into = c(
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
    "gene_revcomp",
    "extra_stuff"
  )) %>%
  select(-accession, -species, -pos, -len, -codonstart, -correction,
         -numAA, -matches("nuc"), -partial, -numchars, -extra_stuff) %>%
  separate(gene, sep = "\\*", into = c("gene", "allele")) %>%
  mutate(regiontype = case_when(regiontype == "EX2R" ~ "EX2",
                                regiontype == "EX2T" ~ "EX2",
                                TRUE ~ regiontype)) %>%
  # don't care about orphons
  filter(!functionality %in% c("ORF")) %>%
  # ignore hinge exons
  filter(!str_starts(regiontype, "H")) %>%
  # ignore membrane and soluble splice sites
  filter(!str_starts(regiontype, "M")) %>%
  filter(regiontype != "CHS") %>%
  # IGHEP isn't consistently documented
  filter(!str_starts(gene, "IGHEP")) %>%
  # not sure what this is (IgG4a?)
  filter(!str_starts(gene, "IGHG4A")) %>%
  # filter out leaders
  filter(!str_starts(gene, "IGLL")) %>%
  # ignore open reading frames
  filter(!str_detect(gene, "/OR")) %>%
  ## append the region site on constant regions (except IGK/L) since these are
  ## documented at the exon level
  mutate(gene = if_else(
    str_detect(gene, "^(IGH(M|D|E.*|G.*|A.*)|TR(A|D|B|G)C.*)$"),
    sprintf("%s-%s", gene, regiontype),
    gene
  )) %>%
  mutate(gene = if_else(str_starts("IGH
  select(-regiontype) %>%
  full_join(order_df, by = "gene")
  ## select(-functionality) %>%
  ## filter(!is.na(allele)) %>%
  ## filter(!is.na(gene_pos))
 
