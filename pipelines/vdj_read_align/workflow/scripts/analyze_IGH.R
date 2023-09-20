library(tidyverse)
library(umap)
library(ggpubr)

df <- read_tsv(
  "../../results/immuno_alignments/ont/split/IGH/final/allele.tsv"
) %>%
  mutate(gene = sprintf("%s%s*%s", major, minor, allele))

subset_major <- function(df, x, usemin = FALSE) {
  .index_col <- paste0(x, "index")
  .pos_col <- paste0(x, "pos")
  .allele_col <- paste0(x, "allele")
  .f <- if (usemin) min else max
  df %>%
    group_by(gene) %>%
    mutate(allele = as.integer(factor(allele))) %>%
    ungroup() %>%
    filter(major == x) %>%
    group_by(rname) %>%
    filter(.f(pos) == pos) %>%
    ungroup() %>%
    select(rname, gene, gene_pos, pos, allele) %>%
    rename(
      !! x := gene,
      !! .index_col := gene_pos,
      !! .pos_col := pos,
      !! .allele_col := allele
    )
}

v_df <- subset_major(df, "V")
dmin_df <- subset_major(df, "D", TRUE) %>%
  rename_with(~ paste0(.x, "_min"), starts_with("D"))
dmax_df <- subset_major(df, "D") %>%
  rename_with(~ paste0(.x, "_max"), starts_with("D"))
j_df <- subset_major(df, "J", TRUE)
c_df <- subset_major(df, "C", TRUE)

all_vdjc_df <- v_df %>%
  full_join(dmin_df, by = "rname") %>%
  full_join(dmax_df, by = "rname") %>%
  full_join(j_df, by = "rname") %>%
  full_join(c_df, by = "rname") %>%
  mutate(
    combo = str_c(V, J, C, sep = "_"),
    Vindex = 159 - Vindex,
    Dlen = Dindex_max - Dindex_min,
    Dindex = Dindex_min - 160,
    Jindex = Jindex - 186,
    Cindex = Cindex - 196,
    ## VDmin_dist = Dpos_min - Vpos,
    ## VDmax_dist = Dpos_max - Vpos,
    ## VJ_dist = Jpos - Vpos,
    ## VC_dist = Cpos - Vpos,
    ## DminDmax_dist = Dpos_max - Dpos_min,
    ## DminJ_dist = Jpos - Dpos_min,
    ## DminC_dist = Cpos - Dpos_min,
    ## DmaxJ_dist = Jpos - Dpos_max,
    ## DmaxC_dist = Cpos - Dpos_max,
    ## JC_dist = Cpos - Jpos
  ) %>%
  select(-Dindex_min, -Dindex_max) %>%
  replace_na(list(Dlen = 0, D_min = "UNK", D_max = "UNK")) %>%
  mutate(across(matches("^(V|D|J|C).+", ignore.case = F), ~ replace_na(.x, -1)))

vdjc_df <- all_vdjc_df %>%
  filter(Vindex > -1 & Cindex > -1)

jc_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  filter(Jindex > -1 & Cindex > -1)

vj_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  filter(Vindex > -1 & Jindex > -1)

vd_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  anti_join(vj_df, by = "rname") %>%
  filter(Vindex > -1 & Dlen > 0)

dj_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  anti_join(vj_df, by = "rname") %>%
  anti_join(vd_df, by = "rname") %>%
  filter(Dlen > 0 & Jindex > -1)

v_only_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  anti_join(vj_df, by = "rname") %>%
  anti_join(vd_df, by = "rname") %>%
  anti_join(dj_df, by = "rname") %>%
  filter(Vindex > -1)

c_only_df <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  anti_join(vj_df, by = "rname") %>%
  anti_join(vd_df, by = "rname") %>%
  anti_join(dj_df, by = "rname") %>%
  anti_join(v_only_df, by = "rname") %>%
  filter(Cindex > -1)

remainder <- all_vdjc_df %>%
  anti_join(vdjc_df, by = "rname") %>%
  anti_join(jc_df, by = "rname") %>%
  anti_join(vj_df, by = "rname") %>%
  anti_join(vd_df, by = "rname") %>%
  anti_join(dj_df, by = "rname") %>%
  anti_join(v_only_df, by = "rname") %>%
  anti_join(c_only_df, by = "rname")

group_list <- function(df) {
  group_map(df, ~ .x) %>% set_names(unlist(group_keys(df)))
}

plot_umap_alleles <- function(df) {
  umap_out <- df %>%
    select(matches("index"), matches("allele"), Dlen) %>%
    select(where(~ any(.x > -1))) %>%
    umap()
  .df <- umap_out$layout %>%
    as_tibble(.name_repair = "minimal") %>%
    set_names(c("V1", "V2")) %>%
    bind_cols(df) %>%
    select(V1, V2, rname, V, D_min, D_max, J, C) %>%
    select(where(~ all(!is.na(.x))))
  .df %>%
    pivot_longer(cols = -c(rname, V1, V2), names_to = "major", values_to = "allele") %>%
    group_by(major) %>%
    group_list() %>%
    imap(~ ggplot(.x, aes(V1, V2, color = allele)) +
           geom_point() +
           labs(x = NULL, y = NULL, color = paste0(.y, " allele"))) %>%
    ggarrange(plotlist = .)
}

plot_umap_alleles(vdjc_df)

plot_umap_alleles(jc_df)

df %>%
  right_join(vd_df, by = "rname") %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()

df %>%
  right_join(vj_df, by = "rname") %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()

df %>%
  right_join(dj_df, by = "rname") %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()

df %>%
  right_join(v_only_df, by = "rname") %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()

df %>%
  right_join(c_only_df, by = "rname") %>%
  ggplot(aes(pos, gene_pos, group = rname, color = major)) +
  geom_point() +
  geom_line()


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

df %>%
  filter(minor == "3-60" | minor == "(III)-38-1") %>%
  pivot_wider(id_cols = rname, values_from = pos, names_from = minor) %>%
  mutate(diff = `(III)-38-1` - `3-60`) %>%
  ggplot(aes(diff)) +
  geom_histogram()
