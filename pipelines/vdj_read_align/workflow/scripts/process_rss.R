library(tidyverse)

flank_df <- readr::read_tsv(
  snakemake@input,
  col_types = "cicidcciicccccc",
  col_names = c(
    "qname",
    "flag",
    "rname",
    "pos",
    "mapq",
    "cigar",
    "rnext",
    "pnext",
    "tlen",
    "seq",
    "qual",
    "XA",
    "MD",
    "NM",
    "XM"
  )) %>%
  mutate(edit = as.numeric(str_replace(NM, "NM:i:", ""))) %>%
  select(qname, flag, rname, pos, edit)

get_rss <- function(df, forward) {
  # filter by forward/reverse
  flag_df <- df %>%
    filter(flag == ifelse(forward, 0, 16)) %>%
    select(-flag)
  # grab all the 5p motifs, and make columns of positions within which the other
  # motif should fall if it is valid
  offset <- ifelse(forward, 7, 9)
  df5p <- flag_df %>%
    filter(qname == ifelse(forward, "heptamer", "nonamer")) %>%
    mutate(lower12 = pos + offset + 10,
           upper12 = pos + offset + 14,
           lower23 = pos + offset + 21,
           upper23 = pos + offset + 25) %>%
    select(-qname)
  # grab all 3p motifs and relabel columns so they don't clash
  df3p <- flag_df %>%
    filter(qname == ifelse(forward, "nonamer", "heptamer")) %>%
    select(-qname) %>%
    rename(pos3p = pos, edit3p = edit)
  # do a fuzzy inner join on the 12bp positions
  df12 <- df5p %>%
    inner_join(df3p, join_by(rname, between(y$pos3p, x$lower12, x$upper12))) %>%
    mutate(is12 = TRUE)
  # ...and same with the 23bp positions
  df23 <- df5p %>%
    inner_join(df3p, join_by(rname, between(y$pos3p, x$lower23, x$upper23))) %>%
    mutate(is12 = FALSE)
  # combine rows and summarize
  bind_rows(df12, df23) %>%
    mutate(edit = edit + edit3p,
           is_reversed = !forward,
           spacer = pos3p - pos - offset,
           end = pos3p + ifelse(forward, 9, 7)) %>%
    select(rname, pos, end, edit, spacer, is12)
}

# write output to a bed file for easy debugging as necessary;
#   name -> indicates spacer size
#   thickStart/End -> indicates spacer start position
#   score -> correlates with edit distance
get_rss(flank_df, TRUE) %>%
  bind_rows(get_rss(flank_df, FALSE)) %>%
  arrange(rname, pos) %>%
  rename("#rname" = rname,
         start = pos) %>%
  mutate(name = if_else(is12, "RSS12", "RSS23"),
         score = 1000 - edit * 200) %>%
  relocate(`#rname`, start, end, name, score, strand, thickStart, thickEnd) %>%
  write_tsv(snakemake@output)
