library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/TRA/final/allele.tsv"
) %>%
  mutate(locus = if_else(gene_pos == 40 | (gene_pos >= 63 & gene_pos <= 72), "TRD", "TRA"),
         gene = sprintf("%s%s%s", locus, major, minor))

v_df <- df %>%
  filter(major == "V") %>%
  filter(locus == "TRA") %>%
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

v_df %>%
  inner_join(d_df, by = "rname") %>%
  inner_join(j_df, by = "rname") %>%
  mutate(combo = str_c(V, D, J, sep = "_")) %>%
  ggplot(aes(y = combo)) +
  geom_bar()
