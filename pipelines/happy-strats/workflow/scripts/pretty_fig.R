library(tidyverse)

root <- "../../results"

pretty_theme <-
  theme(
    text = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.25, "mm"),
    axis.line = element_line(linewidth = 0.15),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0.75, 0, "mm"),
    plot.margin = margin(1, 1, 1, 1, "mm"),
    legend.key.size = unit(0.55, "lines"),
    legend.spacing = unit(1.5, "mm"),
    strip.text = element_text(size = rel(1), margin = margin(1, 1, 1, 1, "mm")),
    strip.background = element_rect(
      linetype = "blank",
      fill = "gray"
    )
  )

geom_line_size <- 0.25
geom_point_size <- 0.4
annot_size <- 6
errorbar_size <- 0.2

phred <- function(x) {
  if_else(x <= 0, NA, - 10 * log10(x))
}

binCI <- function(x, n, name) {
  .mean <- sprintf("%s_mean", name)
  .lower <- sprintf("%s_lower", name)
  .upper <- sprintf("%s_upper", name)
  if_else(
    n == 0, 
    tibble(PointEst = NA, Lower = NA, Upper = NA),
    Hmisc::binconf(x, n, alpha = 0.05, method = "wilson") %>%
      as_tibble()
  ) %>%
    rename(
    {{ .mean }} := PointEst,
    {{ .lower }} := Lower,
    {{ .upper }} := Upper
  )
}

comparison_subsets <- c(
  "AllTandemRepeatsandHomopolymers_slop5",
  "alldifficultregions",
  "alllowmapandsegdupregions",
  "gclt25orgt65_slop50",
  "lowmappabilityall",
  "segdups"
)

not_comparison_subsets <- c(
  "notinalldifficultregions"
)

ont_subsets <- c(
  "AllAutosomes",
  "notinAllHomopolymers_ge7bp_imperfectge11bp_slop5",
  "AllHomopolymers_ge7bp_imperfectge11bp_slop5",
  "notinAllTandemRepeatsandHomopolymers_slop5",
  "AllTandemRepeatsandHomopolymers_slop5"
)

nonsyn_subsets <- c(
  "AllAutosomes",
  "AllTandemRepeatsandHomopolymers_slop5",
  "alldifficultregions",
  "alllowmapandsegdupregions",
  "gclt25orgt65_slop50",
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
    mutate(
      .recall = map2(TRUTH.TP, TRUTH.TP + TRUTH.FN, ~ binCI(.x, .y, "Recall")),
      .precision = map2(QUERY.TP, QUERY.TP + QUERY.FP, ~ binCI(.x, .y, "Precision"))
    ) %>% 
    unnest(c(.recall, .precision)) %>%
    select(Type, Subtype, Subset, matches("^(Recall|Precision)_")) %>%
    mutate(build = which) %>%
    pivot_longer(
      cols = c(-Type, -Subtype, -Subset, -build),
      names_to = c("metric", "bound"),
      names_pattern = "(.+)_(.+)",
      values_to = "value"
    ) %>%
    mutate(
      metric = sprintf("1 - %s", metric),
      value = phred(1 - value)
    ) %>%
    pivot_wider(
      id_cols = c(Type, Subtype, Subset, build, metric),
      names_from = bound,
      values_from = value
    )
}

make_plot <- function(df) {
  df %>%
    mutate(
      Subset = case_when(
        Subset == "AllAutosomes" ~ "Autosomes",
        Subset == "alldifficultregions" ~ "Difficult",
        Subset == "alllowmapandsegdupregions" ~ "Lowmap/Segdup",
        Subset == "gclt25orgt65_slop50" ~ "High/Low GC",
        Subset == "lowmappabilityall" ~ "Lowmap",
        Subset == "segdups" ~ "Segdups",
        Subset == "notinalllowmapandsegdupregions" ~ "Not Lowmap/Segdup",
        Subset == "notinalldifficultregions" ~ "Not Difficult",
        Subset == "notinAllHomopolymers_ge7bp_imperfectge11bp_slop5" ~ "Not HPs",
        Subset == "AllHomopolymers_ge7bp_imperfectge11bp_slop5" ~ "HPs",
        Subset == "notinAllTandemRepeatsandHomopolymers_slop5" ~ "Not TRs/HPs",
        Subset == "AllTandemRepeatsandHomopolymers_slop5" ~ "TRs/HPs",
        TRUE ~ "other"
      )
    ) %>%
    mutate(Type = if_else(Type == "SNP", "SNV", Type)) %>%
    ggplot(aes(mean, fct_rev(Subset), fill = build)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      position = "dodge",
      linewidth = 0.1
    ) +
    theme(legend.position = "top") +
    pretty_theme
}

map_dfr(snakemake@input[["comparison"]], ~ read_summary(path_to_comp, .x)) %>%
  filter(Subset %in% comparison_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(build = factor(build, levels = c("GRCh37", "GRCh38", "CHM13"))) %>%
  make_plot() +
  facet_wrap(c("Type", "metric"), ncol = 2) +
  labs(x = "Phred(Metric)",
       y = NULL,
       fill = "Reference")
ggsave(snakemake@output[["comparison"]], units = "mm", width = 70, height = 90)

map_dfr(snakemake@input[["comparison"]], ~ read_summary(path_to_comp, .x)) %>%
  filter(Subset %in% not_comparison_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(build = factor(build, levels = c("GRCh37", "GRCh38", "CHM13"))) %>%
  make_plot() +
  facet_wrap(c("Type", "metric"), ncol = 2) +
  labs(x = "Phred(Metric)",
       y = NULL,
       fill = "Reference")
ggsave(snakemake@output[["not_comparison"]], units = "mm", width = 70, height = 46)

map_dfr(snakemake@input[["ont"]], ~ read_summary(path_to_ont, .x)) %>%
  filter(Subset %in% ont_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(
    build = if_else(build == "ont_new", "guppy5+clair3", "guppy4+clair1"),
  ) %>%
  make_plot() +
  facet_wrap(c("Type", "metric"), ncol = 2) +
  labs(x = "Phred(Metric)",
       y = NULL,
       fill = "Base/Variant Caller")
ggsave(snakemake@output[["ont"]], units = "mm", width = 70, height = 60)

map_dfr(snakemake@input[["nonsyn"]], ~ read_summary(path_to_ont, .x)) %>%
  filter(Subset %in% nonsyn_subsets) %>%
  filter(Subtype == "*") %>%
  mutate(
    build = if_else(build == "syntenic", "all", "non-syntenic"),
  ) %>%
  make_plot() +
  facet_wrap(c("Type", "metric"), ncol = 2) +
  labs(x = "Phred(Metric)",
       y = NULL,
       fill = "Target Regions")
ggsave(snakemake@output[["nonsyn"]], units = "mm", width = 70, height = 70)
