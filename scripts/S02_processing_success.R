#!/usr/bin/env Rscript
# Supplementary Figure 2: specimen outcomes during COI and 16S processing.
source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
category_order <- c("invalid", "ambiguous", "valid", "inadmissible", "Pre-selection")
category_labels <- c(
  "invalid" = "Contains\ninvalid OTUs",
  "ambiguous" = "Contains only\nambiguous OTUs",
  "valid" = "Valid",
  "inadmissible" = "Rejection\n(Quality filtering)",
  "Pre-selection" = "Pre-selection"
)
category_colors <- c(
  "invalid" = "#CE2929", "ambiguous" = "#FF8B07", "valid" = "#428953",
  "inadmissible" = "#636363", "Pre-selection" = "#5675D6"
)

coi_outcomes <- raw$prey %>%
  count(Status, Category, name = "n") %>%
  add_row(Status = "Contemporary", Category = "Pre-selection", n = 30L) %>%
  add_row(Status = "Museum", Category = "Pre-selection", n = 82L)

# The original processing record supplies these 16S quality-filtering totals;
# they are not derived from the later per-profile read-depth thresholds.
bacterial_outcomes <- tribble(
  ~Status, ~Category, ~n,
  "Contemporary", "Pre-selection", 30L,
  "Contemporary", "inadmissible", 1L,
  "Museum", "Pre-selection", 82L,
  "Museum", "inadmissible", 55L
)
outcomes <- bind_rows(
  mutate(coi_outcomes, marker = "COI"),
  mutate(bacterial_outcomes, marker = "16S")
) %>% mutate(Category = factor(Category, levels = category_order))
stopifnot(nrow(raw$prey) == 112L, !anyNA(outcomes$Category))
write_csv(outcomes, here::here("output", "tables", "S02_processing_outcomes.csv"))

panels <- lapply(c("COI", "16S"), function(marker) {
  ggplot(filter(outcomes, .data$marker == .env$marker), aes(Category, n, fill = Category)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    scale_x_discrete(labels = category_labels) +
    scale_y_continuous(limits = c(0, 90)) +
    scale_fill_manual(values = category_colors) +
    coord_flip() +
    facet_wrap(vars(Status)) +
    labs(x = NULL, y = "Sample size", title = marker) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 13),
          axis.text = element_text(size = 11), strip.text = element_text(size = 11))
})
figure <- ggpubr::ggarrange(panels[[1]], panels[[2]], ncol = 1,
                          labels = c("A", "B"), heights = c(3, 1.8))
save_supplement_figure(figure, "S02_museum_success.png", width = 9, height = 7)
