source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
bacteria <- raw$bacteria
bacteria_taxa <- raw$bacteria_taxa

state <- make_dominance_state(bacteria, bacteria_taxa)
bac_otu <- state$bacteria_otu
cohort <- contemporary_analysis_cohort(state, raw$prey)

# One canonical display map is shared with Figure 4, Figure 5 and S5.
haplotype_sample_map <- make_fig04_arg_map(state$state_tbl, state$cont_ab) %>%
  select(ID, sample_label = arg_label)
stopifnot(
  !anyDuplicated(haplotype_sample_map$ID),
  !anyDuplicated(haplotype_sample_map$sample_label),
  setequal(haplotype_sample_map$ID, state$cont_ab$ID)
)

read_fasta <- function(path) {
  if (!file.exists(path)) {
    return(tibble(OTU = character(), sequence = character()))
  }
  x <- readLines(path, warn = FALSE)
  header_idx <- which(startsWith(x, ">"))
  if (length(header_idx) == 0) {
    return(tibble(OTU = character(), sequence = character()))
  }
  end_idx <- c(header_idx[-1] - 1L, length(x))
  map2_dfr(header_idx, end_idx, function(start, end) {
    header <- str_remove(x[start], "^>")
    seq <- paste0(x[(start + 1L):end], collapse = "")
    tibble(
      OTU = str_to_upper(header),
      sequence = seq
    )
  })
}

fasta <- read_fasta(here::here("data", "raw", "Tim16S.OTu.fas"))

candidate_taxa <- tag_candidate_features(bacteria_taxa) %>%
  mutate(
    symbiont_group = case_when(
      OTU %in% state$otu_dusa ~ "DUSA",
      OTU %in% state$otu_meso ~ "Mesoplasma",
      TRUE ~ "Other candidate"
    )
  ) %>%
  filter(symbiont_group %in% c("DUSA", "Mesoplasma")) %>%
  select(OTU, symbiont_group, tax_text)

feature_totals <- bacteria %>%
  select(ID, Status, Taxon, Year, Country, BacteriaReads, all_of(bac_otu)) %>%
  pivot_longer(cols = all_of(bac_otu), names_to = "OTU", values_to = "reads") %>%
  inner_join(candidate_taxa, by = "OTU")

feature_summary <- feature_totals %>%
  group_by(symbiont_group, OTU, tax_text) %>%
  summarize(
    total_reads_all = sum(reads, na.rm = TRUE),
    n_samples_detected_all = sum(reads > 0, na.rm = TRUE),
    total_reads_historical_abruennichi = sum(
      reads[Status != "Contemporary" & str_detect(str_to_lower(Taxon), "bruennichi")],
      na.rm = TRUE
    ),
    n_historical_abruennichi_detected = sum(
      reads > 0 & Status != "Contemporary" & str_detect(str_to_lower(Taxon), "bruennichi"),
      na.rm = TRUE
    ),
    total_reads_contemporary_abruennichi = sum(
      reads[Status == "Contemporary" & str_detect(str_to_lower(Taxon), "bruennichi")],
      na.rm = TRUE
    ),
    n_contemporary_abruennichi_detected = sum(
      reads > 0 & Status == "Contemporary" & str_detect(str_to_lower(Taxon), "bruennichi"),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  left_join(fasta, by = "OTU") %>%
  mutate(seq_length = nchar(coalesce(sequence, ""))) %>%
  arrange(symbiont_group, desc(total_reads_all), OTU)

readr::write_csv(
  feature_summary,
  here::here("output", "tables", "table_8_symbiont_asv_haplotype_summary.csv")
)

historical_dusa_proportions <- feature_totals %>%
  filter(
    symbiont_group == "DUSA",
    Status != "Contemporary",
    str_detect(str_to_lower(Taxon), "bruennichi")
  ) %>%
  group_by(ID, Taxon, Year, Country, BacteriaReads) %>%
  mutate(
    group_reads = sum(reads, na.rm = TRUE),
    rel_within_dusa = if_else(group_reads > 0, reads / group_reads, NA_real_),
    rel_total_16S = reads / pmax(BacteriaReads, 1)
  ) %>%
  ungroup() %>%
  filter(group_reads > 0) %>%
  left_join(candidate_taxa, by = c("OTU", "symbiont_group", "tax_text")) %>%
  mutate(
    year_num = suppressWarnings(as.numeric(Year)),
    sample_label = paste0(ID, " (", Year, ")")
  ) %>%
  arrange(year_num, ID, desc(reads))

readr::write_csv(
  historical_dusa_proportions,
  here::here("output", "tables", "fig03b_historical_dusa_haplotype_proportions_data.csv")
)

contemporary_symbiont_proportions <- feature_totals %>%
  filter(
    Status == "Contemporary",
    str_detect(str_to_lower(Taxon), "bruennichi")
  ) %>%
  group_by(ID, symbiont_group) %>%
  mutate(
    group_reads = sum(reads, na.rm = TRUE),
    rel_within_group = if_else(group_reads > 0, reads / group_reads, NA_real_)
  ) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(
    focal_symbiont_reads = sum(reads, na.rm = TRUE),
    rel_within_focal_symbionts = if_else(focal_symbiont_reads > 0, reads / focal_symbiont_reads, NA_real_),
    rel_total_16S = reads / pmax(BacteriaReads, 1)
  ) %>%
  ungroup() %>%
  filter(focal_symbiont_reads > 0) %>%
  left_join(state$state_tbl %>% select(ID, dominance_state), by = "ID") %>%
  left_join(haplotype_sample_map, by = "ID") %>%
  left_join(cohort %>% select(ID, S07_included = paired), by = "ID") %>%
  mutate(
    display_group = paste(if_else(symbiont_group == "DUSA", "DUSA-like", symbiont_group), OTU)
  ) %>%
  arrange(symbiont_group, OTU, ID)

readr::write_csv(
  contemporary_symbiont_proportions,
  here::here("output", "tables", "S07_symbiont_haplotypes_data.csv")
)

readr::write_csv(
  tibble(
    metric = c(
      "n_dusa_asv_haplotypes",
      "n_mesoplasma_asv_haplotypes",
      "historical_abruennichi_samples_with_dusa_reads",
      "contemporary_abruennichi_samples_with_focal_symbiont_reads"
    ),
    value = c(
      length(state$otu_dusa),
      length(state$otu_meso),
      n_distinct(historical_dusa_proportions$ID),
      n_distinct(contemporary_symbiont_proportions$ID)
    )
  ),
  here::here("output", "statistics", "09_symbiont_asv_haplotype_summary.csv")
)

if (nrow(contemporary_symbiont_proportions) > 0) {
  sample_order <- state$state_tbl %>%
    filter(ID %in% contemporary_symbiont_proportions$ID) %>%
    left_join(haplotype_sample_map, by = "ID") %>%
    arrange(as.integer(sub("^Arg", "", sample_label))) %>% pull(sample_label)

  readr::write_csv(
    haplotype_sample_map %>% mutate(plot_order = match(sample_label, sample_order)) %>% arrange(plot_order),
    here::here("output", "tables", "S07_symbiont_haplotypes_sample_map.csv")
  )
  # Each bar retains the focal-read denominator used by the original figure.
  stopifnot(
    n_distinct(contemporary_symbiont_proportions$ID) == 30L,
    n_distinct(contemporary_symbiont_proportions$ID[contemporary_symbiont_proportions$S07_included]) == 28L,
    all(abs(contemporary_symbiont_proportions %>%
      group_by(ID) %>% summarize(total = sum(rel_within_focal_symbionts), .groups = "drop") %>%
      pull(total) - 1) < 1e-12)
  )

  plot_data <- contemporary_symbiont_proportions %>%
    mutate(plot_order = match(sample_label, sample_order),
           included = S07_included,
           plotted_fraction = if_else(included, rel_within_focal_symbionts, NA_real_))
  readr::write_csv(plot_data, here::here("output", "tables", "S07_symbiont_haplotypes_plot_data.csv"))

  s08_contemporary_symbiont_asv <- plot_data %>% filter(included) %>%
    mutate(
      sample_label = factor(sample_label, levels = sample_order),
      display_group = factor(
        display_group,
        levels = contemporary_symbiont_proportions %>%
          distinct(symbiont_group, OTU, display_group) %>%
          mutate(otu_number = as.integer(str_remove(OTU, "^OTU"))) %>%
          arrange(symbiont_group, otu_number) %>%
          pull(display_group)
      )
    ) %>%
    ggplot(aes(sample_label, plotted_fraction, fill = display_group)) +
    geom_col(color = "white", linewidth = 0.1, width = 0.85) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = c(
      "DUSA-like OTU2" = "#ff6a0a",
      "Mesoplasma OTU3" = "#eeee00",
      "Mesoplasma OTU145" = "#287A78"
    )) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       breaks = seq(0, 1, .2), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Contemporary A. bruennichi sample",
      y = "Proportion of DUSA-like + Mesoplasma reads",
      fill = "16S sequence variant"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.margin = margin(8, 12, 8, 8)
    )

  # Mark either failed marker in this paired comparison, retaining all30 fixed slots.
  excluded <- cohort %>% filter(!paired) %>% mutate(
    x = match(arg_label, sample_order),
    label = if_else(!eligible_16S, paste0("16S: ", BacteriaReads, " reads"),
                    paste0("Prey: ", PreyReads, " reads")))
  hatch <- bind_rows(lapply(excluded$x, function(x) {
    y0 <- seq(-.85, 1, by = .06); low <- pmax(0, y0); high <- pmin(1, y0 + .85)
    keep <- high > low
    left_edge <- x - .425
    tibble(x = left_edge + low[keep] - y0[keep],
           xend = left_edge + high[keep] - y0[keep], y = low[keep], yend = high[keep])
  }))
  s08_contemporary_symbiont_asv <- s08_contemporary_symbiont_asv +
    geom_rect(data = excluded, aes(xmin = x - .425, xmax = x + .425, ymin = 0, ymax = 1),
      inherit.aes = FALSE, fill = "#f6f6f6", color = "#b0b0b0", linewidth = .3) +
    geom_segment(data = hatch, aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, color = "#b0b0b0", linewidth = .3) +
    geom_label(data = excluded, aes(x = x, y = .5, label = label), inherit.aes = FALSE,
      angle = 90, size = 2.5, fill = "white", color = "#444444", linewidth = 0)
  ggplot2::ggsave(
    here::here("figures", "supplement", "S07_symbiont_haplotypes.png"),
    s08_contemporary_symbiont_asv,
    width = 11,
    height = 5.4,
    dpi = 300,
    bg = "white"
  )
  ggplot2::ggsave(
    here::here("figures", "supplement", "S07_symbiont_haplotypes.pdf"),
    s08_contemporary_symbiont_asv, width = 11, height = 5.4, device = "pdf", useDingbats = FALSE, bg = "white"
  )
}

log_line("S07_symbiont_variants completed")
