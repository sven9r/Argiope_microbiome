# Build main Figure 4 from the bacterial and detailed prey exports by specimen ID.
# Run from the project root: make composition
# This uses the project's R environment and preserves marker-specific cutoffs.
suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
})

merge_args <- commandArgs(trailingOnly = TRUE)
merge_options <- list(root = here::here(), output = NULL, dpi = 300)
while (length(merge_args)) {
  key <- sub("^--", "", merge_args[[1]])
  if (!key %in% names(merge_options) || length(merge_args) < 2) {
    stop("Options: --root PATH --output PATH --dpi NUMBER")
  }
  merge_options[key] <- list(merge_args[[2]])
  merge_args <- merge_args[-c(1, 2)]
}
merge_root <- normalizePath(merge_options$root, mustWork = TRUE)
merge_output <- if (is.null(merge_options$output)) {
  file.path(merge_root, "output", "species_symbionts", "merged_figure4A_prey")
} else merge_options$output
merge_dpi <- as.numeric(merge_options$dpi)
stopifnot(length(merge_dpi) == 1, is.finite(merge_dpi), merge_dpi > 0)
dir.create(merge_output, recursive = TRUE, showWarnings = FALSE)

merge_paths <- c(
  canonical_sample_map = "data/reference/arg_sample_map.csv",
  bacterial_composition = "output/tables/table_3_dominance_state_per_sample.csv",
  detailed_prey_composition = "output/tables/S07_detailed_diet_stack_data.csv",
  prey_sample_set = "output/tables/S07_detailed_diet_sample_set.csv",
  owning_workbook = "data/Raw_data_COI_16S_2.xlsm",
  bacterial_panel_script = "scripts/S05_bacterial_profiles.R",
  prey_panel_script = "scripts/S09_prey_profiles.R",
  source_helpers = "scripts/00_setup.R"
)
if (!all(file.exists(file.path(merge_root, merge_paths)))) {
  stop("Required source exports are missing. Run make bacteria prey first.")
}
merge_read <- function(name) readr::read_csv(file.path(merge_root, merge_paths[[name]]), show_col_types = FALSE)
merge_canonical <- merge_read("canonical_sample_map") %>% select(ID, arg_label)
merge_states <- merge_read("bacterial_composition")
merge_diet <- merge_read("detailed_prey_composition")
merge_prey_set <- merge_read("prey_sample_set")
merge_workbook <- file.path(merge_root, merge_paths[["owning_workbook"]])
merge_bac_raw <- readxl::read_excel(merge_workbook, sheet = "16S_raw_data")
merge_coi_raw <- readxl::read_excel(merge_workbook, sheet = "COI_raw_data")

merge_prey_labels <- merge_diet %>% distinct(ID, arg_label = spider_label)
merge_shared <- inner_join(merge_canonical, merge_prey_labels, by = "ID", suffix = c("_16S", "_prey"))
stopifnot(all(merge_shared$arg_label_16S == merge_shared$arg_label_prey))
merge_labels <- bind_rows(merge_canonical, merge_prey_labels) %>%
  distinct(ID, arg_label) %>%
  arrange(as.integer(sub("^Arg", "", arg_label)))
stopifnot(!anyDuplicated(merge_labels$ID), !anyDuplicated(merge_labels$arg_label))

merge_samples <- merge_labels %>%
  left_join(merge_bac_raw %>% select(ID, Taxon, Status, bacterial_reads = BacteriaReads), by = "ID") %>%
  left_join(merge_coi_raw %>% select(ID, prey_taxon = Taxon, prey_status = Status, prey_reads = PreyReads), by = "ID") %>%
  mutate(
    included_16S = bacterial_reads >= 1500,
    included_prey = prey_reads >= 1693,
    matched = included_16S & included_prey,
    reason_16S_excluded = if_else(included_16S, "", paste0(bacterial_reads, " bacterial reads < 1500")),
    reason_prey_excluded = if_else(included_prey, "", paste0(prey_reads, " prey reads < 1693")),
    original_16S_panel = if_else(included_16S, "Figure 4A", ""),
    original_prey_panel = if_else(included_prey, "Former standalone detailed prey supplement; now Figure 4B", ""),
    union_arg_label = arg_label,
    matched_arg_label = if_else(matched, arg_label, "")
  )
stopifnot(
  all(merge_samples$Taxon == "A. bruennichi"), all(merge_samples$prey_taxon == "A. bruennichi"),
  all(merge_samples$Status == "Contemporary"), all(merge_samples$prey_status == "Contemporary"),
  all(merge_samples$included_16S == (merge_samples$bacterial_reads >= 1500)),
  all(merge_samples$included_prey == (merge_samples$prey_reads >= 1693)),
  setequal(merge_prey_set$ID, merge_diet$ID),
  nrow(merge_samples) == 30, sum(merge_samples$included_16S) == 29,
  sum(merge_samples$included_prey) == 29, sum(merge_samples$matched) == 28
)

# Verify that the existing bacterial exports still agree with their workbook.
merge_otu_columns <- grep("^OTU", names(merge_bac_raw), value = TRUE)
merge_count_check <- merge_bac_raw %>%
  filter(ID %in% merge_samples$ID) %>%
  mutate(
    workbook_total = rowSums(across(all_of(merge_otu_columns)), na.rm = TRUE),
    workbook_dusa = OTU2,
    workbook_meso = OTU3 + OTU145
  ) %>%
  select(ID, workbook_total, workbook_dusa, workbook_meso) %>%
  inner_join(merge_states, by = "ID")
stopifnot(
  nrow(merge_count_check) == 30,
  all(merge_count_check$workbook_total == merge_count_check$total_16S_reads),
  all(merge_count_check$workbook_dusa == merge_count_check$dusa_reads),
  all(merge_count_check$workbook_meso == merge_count_check$mesoplasma_reads)
)
merge_bacterial <- merge_states %>%
  filter(ID %in% merge_samples$ID[merge_samples$included_16S]) %>%
  select(ID, `DUSA-like` = dusa_rel, Mesoplasma = mesoplasma_rel, `Other bacteria` = other_rel) %>%
  pivot_longer(-ID, names_to = "group", values_to = "relative_read_abundance_percent") %>%
  mutate(relative_read_abundance_percent = 100 * relative_read_abundance_percent, marker = "16S")
merge_prey <- merge_diet %>%
  mutate(group = if_else(Order == "Other (<3%)", "Other prey orders", as.character(Order))) %>%
  group_by(ID, group) %>%
  summarize(relative_read_abundance_percent = sum(rel_abund), .groups = "drop") %>%
  mutate(marker = "prey_COI")
merge_bacterial_colors <- c(`Other bacteria` = "#999999", Mesoplasma = "#eeee00", `DUSA-like` = "#ff6a0a")
merge_prey_colors <- c(
  Blattodea = "#6a3d9a", Coleoptera = "#a50013", Diptera = "#fb6449", Hemiptera = "#eeee00",
  Hymenoptera = "#74c476", Lepidoptera = "#6baed6", Orthoptera = "#2878b8", `Other prey orders` = "#999999"
)
stopifnot(all(merge_prey$group %in% names(merge_prey_colors)))
merge_totals <- bind_rows(merge_bacterial, merge_prey) %>%
  group_by(ID, marker) %>% summarize(total = sum(relative_read_abundance_percent), .groups = "drop")
stopifnot(all(abs(merge_totals$total - 100) < 1e-7))

merge_export_samples <- merge_samples %>%
  select(ID, arg_label, Taxon, bacterial_reads, prey_reads, included_16S, included_prey, matched,
         reason_16S_excluded, reason_prey_excluded, original_16S_panel, original_prey_panel, union_arg_label, matched_arg_label)
merge_plot_data <- bind_rows(
  tidyr::crossing(ID = merge_samples$ID, group = names(merge_bacterial_colors)) %>% mutate(marker = "16S"),
  tidyr::crossing(ID = merge_samples$ID, group = names(merge_prey_colors)) %>% mutate(marker = "prey_COI")
) %>%
  left_join(bind_rows(merge_bacterial, merge_prey), by = c("ID", "group", "marker")) %>%
  left_join(merge_export_samples, by = "ID") %>%
  mutate(
    included = if_else(marker == "16S", included_16S, included_prey),
    minimum_reads = if_else(marker == "16S", 1500L, 1693L),
    relative_read_abundance_percent = if_else(included, coalesce(relative_read_abundance_percent, 0), NA_real_),
    marker_order = match(marker, c("16S", "prey_COI")),
    group_order = if_else(marker == "16S", match(group, names(merge_bacterial_colors)), match(group, names(merge_prey_colors)))
  ) %>%
  arrange(match(ID, merge_samples$ID), marker_order, group_order) %>%
  select(ID, arg_label, union_arg_label, matched_arg_label, marker, group, relative_read_abundance_percent, included, minimum_reads)

# Draw excluded profiles as hatched slots, including their read depths.
merge_hatching <- function(excluded, slope = 13, step = 5) {
  bind_rows(lapply(excluded$x, function(x) {
    y0 <- seq(-slope * .73, 100, by = step)
    low <- pmax(0, y0)
    high <- pmin(100, y0 + slope * .73)
    keep <- high > low
    left_edge <- x - .365
    tibble(x = left_edge + (low[keep] - y0[keep]) / slope,
           xend = left_edge + (high[keep] - y0[keep]) / slope, y = low[keep], yend = high[keep])
  }))
}

merge_panel <- function(samples, marker, label_field, show_labels) {
  colors <- if (marker == "16S") merge_bacterial_colors else merge_prey_colors
  levels <- if (marker == "16S") rev(names(colors)) else names(colors)
  positions <- samples %>% mutate(x = row_number())
  data <- merge_plot_data %>% filter(.data$marker == .env$marker, included, ID %in% samples$ID) %>%
    left_join(positions %>% select(ID, x), by = "ID") %>% mutate(group = factor(group, levels = levels))
  excluded <- positions %>%
    filter(if (marker == "16S") !included_16S else !included_prey) %>%
    mutate(count = if (marker == "16S") bacterial_reads else prey_reads,
           label = paste0(format(count, big.mark = ",", scientific = FALSE, trim = TRUE), " reads"))
  panel <- ggplot(data, aes(x, relative_read_abundance_percent, fill = group)) +
    geom_col(width = .73, color = "#333333", linewidth = .15) +
    scale_fill_manual(values = colors, breaks = names(colors),
                      labels = if (marker == "16S") c("Other bacteria", "*Mesoplasma* sp.", "DUSA-like") else names(colors)) +
    scale_x_continuous(breaks = positions$x, labels = positions[[label_field]],
                       limits = c(.3, nrow(samples) + .7), expand = c(0, 0)) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = if (marker == "16S") "A   Bacterial 16S composition" else "B   Prey COI composition",
         x = if (show_labels) "Specimen" else NULL, y = "Relative read abundance (%)",
         fill = if (marker == "16S") "Bacterial group" else "Prey order") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 16, face = "bold", margin = margin(b = 12)),
          axis.text.x = if (show_labels) element_text(size = 8, angle = 0, hjust = .5) else element_blank(),
          axis.ticks.x = element_blank(), axis.text.y = element_text(size = 10),
          axis.title.x = element_text(margin = margin(t = 10)), axis.title.y = element_text(margin = margin(r = 8)),
          legend.position = "right", legend.title = element_text(size = 11),
          legend.text = ggtext::element_markdown(size = 10), legend.key.size = grid::unit(11, "pt"),
          plot.margin = margin(10, 12, 10, 10))
  if (nrow(excluded)) {
    panel <- panel +
      geom_rect(data = excluded, aes(xmin = x - .365, xmax = x + .365, ymin = 0, ymax = 100),
                inherit.aes = FALSE, fill = "#f6f6f6", color = "#b0b0b0", linewidth = .3) +
      geom_segment(data = merge_hatching(excluded), aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, color = "#b0b0b0", linewidth = .3) +
      geom_label(data = excluded, aes(x = x, y = 50, label = label), inherit.aes = FALSE,
                 angle = 90, size = 2.7, fill = "white", color = "#444444", linewidth = 0,
                 label.padding = grid::unit(.08, "lines"))
  }
  panel
}

merge_figure <- function(samples, matched_only = FALSE) {
  label_field <- if (matched_only) "matched_arg_label" else "union_arg_label"
  (merge_panel(samples, "16S", label_field, FALSE) / merge_panel(samples, "prey_COI", label_field, TRUE)) +
    patchwork::plot_layout(heights = c(1, 1), guides = "keep")
}

merge_union_plot <- merge_figure(merge_samples)
merge_matched_plot <- merge_figure(merge_samples %>% filter(matched), matched_only = TRUE)
ggsave(file.path(merge_output, "merged_figure4A_prey_union.png"), merge_union_plot, width = 16, height = 8.6, dpi = merge_dpi, bg = "white")
ggsave(file.path(merge_output, "merged_figure4A_prey_matched.png"), merge_matched_plot, width = 16, height = 8.6, dpi = merge_dpi, bg = "white")
ggsave(file.path(merge_root, "figures/main/fig04_bacteria_prey_composition.pdf"), merge_union_plot,
       width = 16, height = 8.6, device = "pdf", useDingbats = FALSE, bg = "white")
merge_main_dir <- file.path(merge_root, "figures", "main")
dir.create(merge_main_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.copy(file.path(merge_output, "merged_figure4A_prey_union.png"),
                    file.path(merge_main_dir, "fig04_bacteria_prey_composition.png"), overwrite = TRUE))
write_csv(merge_export_samples, file.path(merge_output, "sample_map.csv"), na = "")
write_csv(merge_plot_data, file.path(merge_output, "plot_data.csv"), na = "")

merge_provenance <- list(
  sources = lapply(merge_paths, function(path) list(path = path, sha256 = digest::digest(file = file.path(merge_root, path), algo = "sha256"))),
  script = "scripts/04_bacteria_prey_composition.R", cutoffs = list(`16S` = 1500, prey_COI = 1693),
  sample_counts = list(union = nrow(merge_samples), `16S` = sum(merge_samples$included_16S), prey_COI = sum(merge_samples$included_prey), matched = sum(merge_samples$matched)),
  validation = list(canonical_labels_unique = TRUE, shared_labels_match_by_ID = TRUE,
                    relative_abundances_sum_to_100 = TRUE, inclusion_matches_raw_workbook_cutoffs = TRUE,
                    exported_16S_counts_match_workbook = TRUE),
  missing_marker_inclusion = merge_export_samples %>% filter(!matched) %>%
    select(ID, arg_label, bacterial_reads, prey_reads, reason_16S_excluded, reason_prey_excluded),
  software = list(R = R.version.string, ggplot2 = as.character(packageVersion("ggplot2")), patchwork = as.character(packageVersion("patchwork"))),
  notes = c("The detailed prey panel is now Figure 4B; Mantel is Supplementary Figure 9. Legacy S07 data filenames are retained for provenance.",
            "Panels use per-specimen read percentages from the existing Figure 4A and S07 R exports; no averaging across specimens.",
            "Specimens are joined by ID and drawn in the same order in both panels.",
            "Display labels follow focal DUSA fraction and remain fixed in the matched companion; current and previous identities are retained in data/reference/arg_sample_map.csv.",
            "Residual categories are Other bacteria and Other prey orders; canonical colors and Blattodea are retained.")
)
jsonlite::write_json(merge_provenance, file.path(merge_output, "provenance.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("Figure 4: 30 specimen positions, 29 profiles per marker, 28 matched.\n")
