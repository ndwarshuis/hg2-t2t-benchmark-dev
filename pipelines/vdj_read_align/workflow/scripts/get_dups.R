suppressMessages(library(tidyverse))

order_df <- read_tsv(snakemake@input[["order"]], col_types = "ci")

genes_with_rank <- order_df %>% pull(gene)

df <- readr::read_tsv(snakemake@input[["distances"]],
                      col_types = cols(gene = "c", .default = "i"))

dist_df <- df %>%
  pivot_longer(cols = -gene, values_to = "distance", names_to = "other") %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = c(gene, other), values_to = "name", names_to = "var") %>%
  mutate(var = if_else(var == "gene", "geneA", "geneB")) %>%
  mutate(major = str_sub(name, 4, 4),
         minor = str_sub(str_extract(name, "(.*)\\*.*", 1), 5),
         allele = as.integer(str_extract(name, ".*\\*(.*)", 1))) %>%
  filter(str_extract(name, "(.*)\\*.*", 1) %in% genes_with_rank) %>%
  pivot_wider(id_cols = c(id, distance),
              names_from = var,
              values_from = c(name, major, minor, allele),
              names_glue = "{var}_{.value}") %>%
  filter(!is.na(geneA_name)) %>%
  filter(!is.na(geneB_name)) %>%
  select(-id)

dist_df %>%
  filter(geneA_minor != geneB_minor) %>%
  filter(distance < 10) %>%
  ggplot(aes(geneA_name, geneB_name, fill = distance)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(snakemake@output[["plot"]])

sim_df <- dist_df %>%
  filter(geneA_minor != geneB_minor) %>%
  filter(distance < 10) %>%
  arrange(geneA_name) %>%
  mutate(rankA = as.numeric(factor(geneA_name))) %>%
  arrange(geneB_name) %>%
  mutate(rankB = as.numeric(factor(geneB_name))) %>%
  filter(rankA < rankB) %>%
  select(geneA_name, geneB_name, rankA, rankB)

similar_genes <- sim_df %>%
  anti_join(select(sim_df, rankB) %>% unique(), by = c("rankA" = "rankB")) %>%
  group_by(geneA_name) %>%
  group_map(~ c(.y$geneA_name, .x$geneB_name))

partition <- function(xs, fun) {
  # gross...
  list(pass = keep(xs, fun), fail = discard(xs, fun))
}

pick_next <- function(xs, fun, more_funs = list()) {
  pick_first <- function(xs) {
    ys <- sort(xs)
    list(rep = ys[[1]], drop = ys[-1])
  }
  next_fun <- if (length(more_funs) == 0) {
    pick_first
  } else {
    function(ys) pick_next(ys, more_funs[[1]], more_funs[-1])
  }
  res <- partition(xs, fun)
  if (length(res$pass) == 1) {
    list(rep = res$pass[[1]], drop = res$fail)
  } else if (length(res$pass) == 0) {
    next_fun(xs)
  } else {
    .res <- next_fun(res$pass)
    .res$drop = c(.res$drop, res$fail)
    .res
  }
}

is_not_dup <- function(gene) {
  str_sub(str_extract(gene, "(.*)\\*.*", 1), -1) != "D"
}

is_not_K_dup <- function(gene) {
  !str_detect(str_sub(gene, 4), "D-")
}

similar_genes %>%
  map(~ pick_next(.x, ~ !str_detect(.x, "OR"), list(is_not_dup, is_not_K_dup))) %>%
  imap_dfr(~ tibble(id = .y, new = .x$rep, old = .x$drop)) %>%
  bind_rows(tibble(id = integer(), new = character(), old = character())) %>%
  write_tsv(snakemake@output[["tsv"]])
