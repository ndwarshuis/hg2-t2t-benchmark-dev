library(tidyverse)

order_df <- read_tsv(
  "../../static/gene_order.tsv",
  col_types = "ci"
) %>%
  mutate(major = str_sub(gene, 4, 4)) %>%
  filter(major %in% c("V", "D", "J")) %>%
  filter(gene != "IGHD") %>%
  select(-major)

cov_df <- read_tsv(
  "../../results/vdj_mapping/test1/coverage.txt",
  col_types = "ciiiidddd"
) %>%
  mutate(gene = str_extract(`#rname`, ".+\\|(.*)\\|.+", 1)) %>%
  select(-`#rname`) %>%
  separate(col = gene, into = c("gene", "allele"), sep = "\\*") %>%
  mutate(gene = str_replace_all(gene, "[()]", "")) %>%
  left_join(order_df, by = "gene") %>%
  mutate(locus = str_sub(gene, 1, 3),
         major = str_sub(gene, 4, 4),
         minor = str_sub(gene, 5),
         locus = if_else(locus == "TRD", "TRA", locus))

cov_df %>%
  filter(numreads > 0) %>%
  group_by(pos, locus, major, minor) %>%
  summarize(totreads = sum(numreads),
            meandepth = meandepth * numreads / totreads,
            .groups = "drop") %>%
  ggplot(aes(pos, meandepth, fill = major)) +
  geom_col() +
  facet_wrap(c("locus"))
  
           
