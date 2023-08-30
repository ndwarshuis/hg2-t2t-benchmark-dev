library(tidyverse)

df <- read_tsv(
  "../../results/immuno_alignments/hifi/split/IGK/final/allele.tsv"
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

df %>%
  filter(nC > 0) %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_line() +
  geom_point() +
  labs(x = "read position", y = "gene rank", color = "region type")

vjc_df %>%
  ggplot(aes(y = combo)) +
  geom_bar()
# hmm, one without J4, probably just noise (igv confirms this)
# hap1 = V3-20_J4_CDEL

# what happened to the other clone? (plenty of reads without J)

vc_df <- v_df %>%
  anti_join(vjc_df) %>%
  inner_join(c_df, by = "rname") %>%
  mutate(combo = str_c(V, C, sep = "_"))

vc_df %>%
  ggplot(aes(y = combo)) +
  geom_bar()
# again, the one clone that didn't map is just noise (igv confirmed)
# hap2 = V2-29_CDEL

# what about all these reads with Vs that end at gene position 70?
df %>%
  filter(nC == 0, nJ == 0) %>%
  pull(minor) %>%
  unique() %>%
  sort()
# my guess is that these are all from the duplicate end of the locus, since
# each of these genes has a mirror analog, this also agrees with the fact that
# gene position 70 is 3-7, and 3D-7 is the last gene on the mirror half
  
