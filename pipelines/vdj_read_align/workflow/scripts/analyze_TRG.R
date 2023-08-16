library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/TRG/final/allele.tsv"
) %>%
  mutate(gene = sprintf("%s%s", major, minor))

# plot weird deletion...
df %>%
  filter(gene_pos == 7 | gene_pos == 2) %>%
  pivot_wider(id_cols = rname, names_from = gene, values_from = pos) %>%
  mutate(distance = V5P - V2) %>%
  ggplot(aes(distance)) +
  geom_histogram()
