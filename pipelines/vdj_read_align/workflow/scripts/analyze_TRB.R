library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/TRB/final/allele.tsv"
) %>%
  mutate(gene = sprintf("%s%s", major, minor))

v_df <- df %>%
  filter(major == "V") %>%
  group_by(rname) %>%
  filter(max(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(V = gene)

d_df <- df %>%
  filter(major == "D") %>%
  group_by(rname) %>%
  filter(max(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(D = gene)

j_df <- df %>%
  filter(major == "J") %>%
  group_by(rname) %>%
  filter(min(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(J = gene)

vdj_df <- inner_join(v_df, j_df, by = "rname") %>%
  inner_join(d_df, by = "rname")
