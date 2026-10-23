#!/usr/bin/env Rscript
# Maintain the reviewed S2 design from script_analysis_thesis.Rmd:310-369.
# Only panel B's rejection label changes; all bar data and positions are retained.
# Run with the project's installed R 4.5 environment, from the project root.
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
})

root <- normalizePath(getwd())
workbook <- file.path(root, "data/Raw_data_COI_16S_2.xlsm")
target <- file.path(root, "figures/supplement/S02_museum_success.png")
support <- file.path(root, "output/tables")
dir.create(support, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
prey <- read_excel(workbook, sheet = "COI_raw_data")
bacteria <- read_excel(workbook, sheet = "16S_raw_data")
stopifnot(nrow(prey) == 112L, nrow(bacteria) == 112L,
          !anyDuplicated(prey$ID), !anyDuplicated(bacteria$ID))

prey_counts <- prey %>% select(ID, Status, Category) %>%
  group_by(Status, Category) %>% summarize(n = n(), .groups = "drop") %>%
  add_row(Status = "Contemporary", Category = "Pre-selection", n = 30) %>%
  add_row(Status = "Museum", Category = "Pre-selection", n = 82)
expected_prey <- c("Contemporary/inadmissible" = 1, "Contemporary/valid" = 29,
  "Contemporary/Pre-selection" = 30, "Museum/inadmissible" = 52,
  "Museum/valid" = 19, "Museum/ambiguous" = 8, "Museum/invalid" = 3,
  "Museum/Pre-selection" = 82)
stopifnot(nrow(prey_counts) == length(expected_prey),
  all(prey_counts$n == unname(expected_prey[paste(prey_counts$Status, prey_counts$Category, sep = "/")])))

success_labels <- c("Contains\ninvalid OTUs", "Contains only\nambiguous OTUs",
  "Valid", "Rejection\n(Quality filtering)", "Pre-selection")
prey_success <- ggplot(prey_counts) +
  aes(x = reorder(Category, n), y = n, fill = Category) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  geom_col(colour = "black", linewidth = 0.2) + ylim(0, 90) +
  labs(x = " ", y = "\nSample size", title = "COI") + coord_flip() +
  theme_minimal() + theme(legend.position = "none",
    plot.title = element_text(size = 13), plot.caption = element_text(size = 12),
    axis.text = element_text(size = 11), strip.text.x = element_text(size = 11)) +
  facet_wrap(vars(Status)) + scale_x_discrete(labels = success_labels) +
  scale_fill_manual(values = c("#FF8B07", "#636363", "#CE2929", "#5675D6", "#428953"))

# Retain the reviewed source counts, while verifying their threshold correspondence.
bac_counts <- tibble(status = c("Contemporary", "Contemporary", "Museum", "Museum"),
  category = rep(c("Pre-selection", "Rejection\n(Quality filtering)"), 2),
  n = c(30, 1, 82, 55))
stopifnot(sum(bacteria$Status == "Contemporary" & bacteria$BacteriaReads < 507) == 1L,
  sum(bacteria$Status == "Museum" & bacteria$BacteriaReads < 507) == 55L)
bac_success <- ggplot(bac_counts,
  aes(x = reorder(category, n), y = n, fill = category)) +
  geom_col(width = 0.75, colour = "black", linewidth = 0.2) +
  geom_text(aes(label = n), hjust = -0.2, size = 3) + ylim(0, 90) +
  labs(x = " ", y = "\nSample size", title = "16S") + coord_flip() +
  theme_minimal() + theme(legend.position = "none",
    plot.title = element_text(size = 13), plot.caption = element_text(size = 12),
    axis.text = element_text(size = 11), strip.text.x = element_text(size = 11)) +
  facet_wrap(vars(status)) + scale_fill_manual(values = c("#5675D6", "#636363"))

# Edit the text after the original layout has been measured. This preserves every
# gtable width/height, bar coordinate and count-label position, including panel A.
old_label <- "Rejection\n(Quality filtering)"
new_label <- "Below 507\nbacterial reads"
replacements <- 0L
replace_bacterial_axis_label <- function(g) {
  if (inherits(g, "text") && !is.null(g$label)) {
    hit <- !is.na(g$label) & g$label == old_label
    replacements <<- replacements + sum(hit)
    g$label[hit] <- new_label
  }
  if (!is.null(g$grobs)) g$grobs <- lapply(g$grobs, replace_bacterial_axis_label)
  if (!is.null(g$children)) {
    for (i in seq_along(g$children))
      g$children[[i]] <- replace_bacterial_axis_label(g$children[[i]])
  }
  g
}
original_b <- ggplotGrob(bac_success)
labelled_b <- replace_bacterial_axis_label(original_b)
stopifnot(replacements == 1L, identical(original_b$widths, labelled_b$widths),
          identical(original_b$heights, labelled_b$heights),
          identical(original_b$layout, labelled_b$layout))
plot_s02 <- function(b) ggarrange(prey_success, b, ncol = 1,
  labels = c("A", "B"), heights = c(3, 1.8))
reference_path <- tempfile("S02_original_label_", fileext = ".png",
  tmpdir = Sys.getenv("TMPDIR", unset = "/tmp"))
ggsave(reference_path, plot_s02(original_b), width = 9, height = 7,
       dpi = 300, bg = "white")
ggsave(target, plot_s02(labelled_b), width = 9, height = 7,
       dpi = 300, bg = "white")

plot_data <- bind_rows(prey_counts %>% transmute(panel = "A", marker = "COI",
  status = Status, source_category = Category, n),
  bac_counts %>% transmute(panel = "B", marker = "16S", status,
    source_category = category, n))
readr::write_csv(plot_data, file.path(support, "S02_processing_outcome_plot_data.csv"))
readr::write_csv(plot_data %>% transmute(Status = status,
  Category = if_else(marker == "16S" & source_category != "Pre-selection",
                     "inadmissible", source_category), n, marker),
  file.path(support, "S02_processing_outcomes.csv"))
audit <- list(figure = "figures/supplement/S02_museum_success.png",
  change = "Panel B axis label only: Below 507 bacterial reads",
  original_label = old_label, replacement_label = new_label,
  text_grobs_replaced = replacements, width_pixels = 2700L, height_pixels = 2100L,
  original_layout_preserved = TRUE, COI_panel_unchanged = TRUE,
  bacterial_counts = as.list(bac_counts$n),
  original_label_reference_render = reference_path,
  source = "scripts/_archive_20260601/script_analysis_thesis.Rmd:310-369",
  workbook_sha256 = digest::digest(file = workbook, algo = "sha256"),
  output_sha256 = digest::digest(file = target, algo = "sha256"))
jsonlite::write_json(audit, file.path(support, "S02_label_update_audit.json"),
  pretty = TRUE, auto_unbox = TRUE)
cat("S2 updated:", target, "\nOriginal-label reference:", reference_path, "\n")
