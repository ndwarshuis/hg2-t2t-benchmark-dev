library(tidyverse)

root <- "../../results"

subsets <- c(
  "AllAutosomes",
  "AllTandemRepeatsandHomopolymers_slop5",
  "alldifficultregions",
  "alllowmapandsegdupregions",
  "gclt25orgt65_slop50",
  "gclt30orgt55_slop50",
  "lowmappabilityall",
  "segdups"
)

read_summary <- function(path) {
  which = path %>% dirname() %>% dirname() %>% basename()
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
        Subtype %in% c("*", "I16_PLUS", "D16_PLUS") &
        Subset %in% subsets
    ) %>%
    rename(
      Recall = METRIC.Recall,
      Precision = METRIC.Precision,
      ## F1 = METRIC.F1_Score,
      ## Frac_NA = METRIC.Frac_NA
    ) %>%
    select(Type, Subtype, Subset, Recall, Precision) %>%
    mutate(build = which) %>%
    pivot_longer(cols = c(-Type, -Subtype, -Subset, -build),
                 names_to = "metric",
                 values_to = "value") %>%
    mutate(value = - 10 * log10(1 - value))
    ## mutate(value = 1 - value)
}

df <- c(
  "/home/njd2/Dvl/giab/hg2/hg2-t2t-benchmark-dev/pipelines/happy-strats/results/CHM13/comparison/happy/happy.extended.csv",
  "/home/njd2/Dvl/giab/hg2/hg2-t2t-benchmark-dev/pipelines/happy-strats/results/CHM13/non_syntenic/happy/happy.extended.csv"
) %>%
  map_dfr(read_summary)

df %>%
  filter(Subtype == "*") %>%
  ## mutate(build = if_else(build == "ont_new", "guppy4+clair1", "guppy5+clair3")) %>%
  ggplot(aes(value, fct_rev(Subset), fill = build)) +
  geom_col(position = "dodge") +
  facet_wrap(c("Type", "metric"), ncol = 4) +
  labs(x = "Phred(Error Rate)",
       y = NULL,
       fill = "Reference") +
  theme(legend.position = "top") +
  scale_x_continuous(labels = scales::percent)

## ggsave(snakemake@output[[1]], units = "cm", width = 20, height = 8)
