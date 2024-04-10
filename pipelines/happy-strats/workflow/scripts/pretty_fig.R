library(tidyverse)

# TODO use binconf

root <- "../../results"

comparison_subsets <- c(
  ## "AllAutosomes",
  "AllTandemRepeatsandHomopolymers_slop5",
  "alldifficultregions",
  "alllowmapandsegdupregions",
  "gclt25orgt65_slop50",
  "gclt30orgt55_slop50",
  "lowmappabilityall",
  "segdups",
  "notinalllowmapandsegdupregions",
  "notinalldifficultregions"
)

# add (notin)homopolymers and notinalltandemrepeatsandhomopoolymers since that's what ONT is bad at
ont_subsets <- c(
  "AllAutosomes",
  "AllTandemRepeatsandHomopolymers_slop5",
  "alldifficultregions",
  "alllowmapandsegdupregions",
  "gclt25orgt65_slop50",
  "gclt30orgt55_slop50",
  "lowmappabilityall",
  "segdups"
)

path_to_comp <- . %>% dirname() %>% dirname() %>% dirname() %>% basename()
path_to_ont <- . %>% dirname() %>% dirname() %>% basename()

read_summary <- function(f, path) {
  which <- f(path)
  read_csv(
    path,
    col_types = cols(
      Type = "c",
      Subtype = "c",
      Subset = "c",
      Filter = "c",
      Genotype = "c",
      QQ.Field = "c",
      QQ = "c",
      .default = "d"
    ),
    na = c(".", "")
  ) %>%
    filter(
      Filter == "PASS" &
        Subtype %in% c("*", "I16_PLUS", "D16_PLUS")
    ) %>%
    rename(
      Recall = METRIC.Recall,
      Precision = METRIC.Precision,
    ) %>%
    select(Type, Subtype, Subset, Recall, Precision) %>%
    mutate(build = which) %>%
    pivot_longer(cols = c(-Type, -Subtype, -Subset, -build),
                 names_to = "metric",
                 values_to = "value")
}

map_dfr(snakemake@input[["comparison"]], ~ read_summary(path_to_comp, .x)) %>%
  filter(Subset %in% comparison_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(value = 1 - value) %>%
  ggplot(aes(value, fct_rev(Subset), fill = build)) +
  geom_col(position = "dodge") +
  facet_wrap(c("Type", "metric"), scales = "free_x", ncol = 4) +
  labs(x = "Error Rate",
       y = NULL,
       fill = "Reference") +
  theme(legend.position = "top") +
  scale_x_continuous(labels = scales::percent)
ggsave(snakemake@output[["comparison"]], units = "cm", width = 20, height = 8)

map_dfr(snakemake@input[["ont"]], ~ read_summary(path_to_ont, .x)) %>%
  filter(Subset %in% ont_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(value = - 10 * log10(1 - value)) %>%
  ggplot(aes(value, fct_rev(Subset), fill = build)) +
  geom_col(position = "dodge") +
  facet_wrap(c("Type", "metric"), scales = "free_x", ncol = 4) +
  labs(x = "Phred(Error Rate)",
       y = NULL,
       fill = "Base/Variant Caller") +
  theme(legend.position = "top")
ggsave(snakemake@output[["ont"]], units = "cm", width = 20, height = 8)
