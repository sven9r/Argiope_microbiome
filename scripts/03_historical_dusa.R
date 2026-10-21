source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
bacteria <- raw$bacteria
bacteria_taxa <- raw$bacteria_taxa

state <- make_dominance_state(bacteria, bacteria_taxa)
bac_otu <- state$bacteria_otu

summarize_historical_dusa <- function(data, sample_group) {
  if (nrow(data) == 0) {
    return(tibble(
      ID = character(),
      Taxon = character(),
      Year = character(),
      Country = character(),
      Location = character(),
      BacteriaReads = double(),
      DUSA_reads = double(),
      Mesoplasma_reads = double(),
      Other_reads = double(),
      DUSA_rel = double(),
      Mesoplasma_rel = double(),
      sample_group = character()
    ))
  }

  data %>%
    mutate(Country = if_else(Taxon == "A. sector", "Senegal", Country)) %>%
    select(ID, Taxon, Year, Country, any_of("Location"), BacteriaReads, all_of(bac_otu)) %>%
    pivot_longer(cols = all_of(bac_otu), names_to = "OTU", values_to = "reads") %>%
    group_by(ID, Taxon, Year, Country, Location, BacteriaReads) %>%
    summarize(
      DUSA_reads = sum(reads[OTU %in% state$otu_dusa], na.rm = TRUE),
      Mesoplasma_reads = sum(reads[OTU %in% state$otu_meso], na.rm = TRUE),
      Other_reads = sum(reads, na.rm = TRUE) - DUSA_reads - Mesoplasma_reads,
      .groups = "drop"
    ) %>%
    mutate(
      DUSA_rel = DUSA_reads / pmax(BacteriaReads, 1),
      Mesoplasma_rel = Mesoplasma_reads / pmax(BacteriaReads, 1),
      sample_group = sample_group
    )
}

country_code <- function(country) {
  recode(
    country,
    Russia = "ru",
    Greece = "gr",
    Germany = "de",
    France = "fr",
    Senegal = "sn",
    .default = NA_character_
  )
}

country_flag <- function(country) {
  recode(
    country,
    Russia = "\U0001F1F7\U0001F1FA",
    Greece = "\U0001F1EC\U0001F1F7",
    Germany = "\U0001F1E9\U0001F1EA",
    France = "\U0001F1EB\U0001F1F7",
    Senegal = "\U0001F1F8\U0001F1F3",
    .default = country
  )
}

museum_ab <- bacteria %>%
  filter(Status != "Contemporary", str_detect(str_to_lower(Taxon), "bruennichi"))
other_museum_argiope <- bacteria %>%
  filter(
    Status != "Contemporary",
    !str_detect(str_to_lower(Taxon), "bruennichi"),
    str_detect(str_to_lower(Taxon), "^(argiope|a\\.)")
  )

hist_ab_tbl <- summarize_historical_dusa(museum_ab, "museum_A_bruennichi")
hist_other_tbl <- summarize_historical_dusa(other_museum_argiope, "other_museum_Argiope")

readr::write_csv(hist_ab_tbl, here::here("output", "tables", "table_4_historical_dusa_abruennichi.csv"))
readr::write_csv(hist_other_tbl, here::here("output", "tables", "supp_table_6_non_bruennichi_dusa_detection.csv"))

readr::write_csv(
  tibble(
    metric = c(
      "museum_A_bruennichi_samples",
      "museum_A_bruennichi_with_DUSA_reads",
      "other_museum_Argiope_samples",
      "other_museum_Argiope_with_DUSA_reads"
    ),
    value = c(
      nrow(hist_ab_tbl),
      sum(hist_ab_tbl$DUSA_reads > 0, na.rm = TRUE),
      nrow(hist_other_tbl),
      sum(hist_other_tbl$DUSA_reads > 0, na.rm = TRUE)
    )
  ),
  here::here("output", "statistics", "08_historical_dusa_summary.csv")
)

historical_dusa_timeline <- bind_rows(hist_ab_tbl, hist_other_tbl) %>%
  filter(
    DUSA_reads > 0,
    str_detect(str_to_lower(Taxon), "bruennichi") | ID == "TH73"
  ) %>%
  mutate(
    year_num = suppressWarnings(as.numeric(Year)),
    year_label = as.character(as.integer(year_num)),
    flag_code = country_code(Country),
    flag = country_flag(Country),
    sector_label = if_else(
      ID == "TH73",
      "italic('A. sector')",
      NA_character_
    )
  ) %>%
  filter(!is.na(year_num)) %>%
  arrange(year_num, ID) %>%
  mutate(year_label = factor(year_label, levels = unique(year_label)))

readr::write_csv(
  historical_dusa_timeline,
  here::here("output", "tables", "fig03_historical_dusa_timeline_data.csv")
)

if (!requireNamespace("ggflags", quietly = TRUE)) {
  stop("Package 'ggflags' is required for fig03_historical_dusa_timeline.png")
}

s07_historical_dusa <- ggplot(
  historical_dusa_timeline,
  aes(year_label, DUSA_rel)
) +
  ggflags::geom_flag(
    aes(country = flag_code, size = BacteriaReads),
    show.legend = FALSE
  ) +
  geom_point(
    aes(size = BacteriaReads),
    alpha = 0,
    show.legend = TRUE
  ) +
  geom_text_repel(
    aes(label = sector_label),
    parse = TRUE,
    min.segment.length = 0,
    na.rm = TRUE,
    size = 3.2,
    box.padding = 0.25,
    point.padding = 0.6,
    nudge_x = 0.32,
    nudge_y = 0.07,
    seed = 12,
    show.legend = FALSE
  ) +
  scale_size_continuous(range = c(3, 9), breaks = c(100, 1000, 2000, 3000)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "Collection year",
    y = "DUSA-like relative abundance",
    size = "16S reads"
  ) +
  guides(
    size = guide_legend(
      override.aes = list(alpha = 1, color = "black", fill = "black", shape = 16)
    )
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

save_main_figure(
  s07_historical_dusa,
  "fig03_historical_dusa_timeline.png",
  width = 8,
  height = 5
)

log_line("03_historical_dusa completed")
