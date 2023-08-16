library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/IGL/final/allele.tsv"
) %>%
  mutate(gene = sprintf("%s%s", major, minor))

v_df <- df %>%
  filter(major == "V") %>%
  group_by(rname) %>%
  filter(max(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(V = gene)

c_df <- df %>%
  filter(major == "C") %>%
  group_by(rname) %>%
  filter(min(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(C = gene)

vjc_df <- df %>%
  filter(major == "J") %>%
  group_by(rname) %>%
  filter(min(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(J = gene) %>%
  inner_join(v_df, by = "rname") %>%
  inner_join(c_df, by = "rname") %>%
  mutate(combo = str_c(V, J, C, sep = "_"))

vjc_df %>%
  ggplot(aes(y = combo)) +
  geom_bar()
# unarranged hap: V3-1_J1_C1
# arranged hap: V2-14_J3_C3
#
# the rest is just noise
# J6/C3 is impossible (in this case J3 didn't map)
# all other cases are genes that mapped immediately upstream of the dominate hap
# genes (most of these are pseudogenes so couldn't even happen for the rearranged
# hap)
