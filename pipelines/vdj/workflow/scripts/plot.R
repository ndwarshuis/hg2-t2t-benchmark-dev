library(tidyverse)

df_loci <- readr::read_tsv(
  snakemake@input$loci,
  col_types = cols(
    gene = "c",
    global_start = "i",
    global_end = "i",
    hap = "c",
  )
)

df_genes <- readr::read_tsv(
  snakemake@input$genes,
  col_types = cols(
    pos = "i",
    mapped_pos = "i",
    correct_map = "l",
    hap = "c",
    locus = "c"
  )
)

df_CHM13 <- df_genes %>%
  filter(!is.na(pos)) %>%
  filter(!is.na(correct_map)) %>%
  group_by(locus, hap) %>%
  arrange(hap, locus, pos) %>%
  mutate(pos_CHM13 = row_number()) %>%
  ungroup()

df_CHM13 %>%
  ggplot(aes(pos_CHM13, pos_real, color = correct_map, shape = major)) +
  geom_point() +
  facet_wrap(c("locus", "hap"), scales = "free")
ggsave(snakemake@output$unmapped, width = 10, height = 7)

## df_CHM13 %>%
##   filter(locus %in% c("TRA", "TRD")) %>%
##   ggplot(aes(pos_CHM13, pos_real, color = correct_map, shape = major)) +
##   geom_point() +
##   facet_wrap(c("locus", "hap"), scales = "free")

## df_CHM13 %>%
##   filter(locus %in% c("TRB")) %>%
##   ggplot(aes(pos_CHM13, pos_real, color = correct_map, shape = major)) +
##   geom_point() +
##   facet_wrap(c("locus", "hap"), scales = "free")

# genes that are mapped

df_mapped <- df_CHM13 %>%
  filter(correct_map) %>%
  left_join(df_loci, by = c("locus" = "gene", "hap" = "hap")) %>%
  mutate(global_contains = global_start <= start & end <= global_end)

df_landmarks <- df_CHM13 %>%
  filter(locus == "IGK", major == "V", minor == "2-40") %>%
  mutate(landmark = "hsat2")

df_mapped %>%
  ## filter(!is.na(global_contains)) %>%
  ggplot(aes(mapped_pos, pos_CHM13, color = major)) +
  geom_hline(data = df_landmarks, aes(yintercept = pos_CHM13)) +
  geom_point() +
  facet_wrap(c("locus", "hap"), scales = "free")
ggsave(snakemake@output$mapped, width = 10, height = 7)
 
