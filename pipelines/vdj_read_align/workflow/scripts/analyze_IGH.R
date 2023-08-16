library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/IGH/final/allele.tsv"
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

vjc_df <- df %>%
  filter(major == "C") %>%
  group_by(rname) %>%
  filter(min(gene_pos) == gene_pos) %>%
  select(rname, gene) %>%
  rename(C = gene) %>%
  inner_join(v_df, by = "rname") %>%
  inner_join(j_df, by = "rname") %>%
  mutate(combo = str_c(V, J, C, sep = "_"))

vjc_df %>%
  ggplot(aes(y = combo)) +
  geom_bar()

# investigate the non-M clone(s) further
nonM_df <- vjc_df %>%
  filter(C != "CM") %>%
  select(rname) %>%
  left_join(df, by = "rname") %>%
  filter(gene %in% c("CG3", "CG1", "J3", "V3-48", "V3-49", "V(II)-49-1")) %>%
  select(rname, gene, pos)

# plotting distances shows that the weirdness in the bar chart is due to mismapping
full_join(nonM_df, nonM_df, by = "rname", relationship = "many-to-many") %>%
  filter(gene.x != gene.y) %>%
  mutate(distance = pos.y - pos.x) %>%
  mutate(key = map2_chr(gene.x, gene.y, ~ do.call(paste, as.list(sort(c(.x, .y)))))) %>%
  group_by(rname, key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  ggplot(aes(distance, key)) +
  geom_jitter()
# only one recombination: V3-48_J3_CG3 (with some unknown D that can't possibly map)

# investigate the M clone(s) further for weirdness in D
vjc_df %>%
  filter(C == "CM") %>%
  inner_join(d_df, by = "rname") %>%
  mutate(vdjc_combo = str_c(V, D, J, C, sep = "_")) %>%
  group_by(vdjc_combo) %>%
  tally()
# only one recombination: V6-1_D3-22_J3_CM
