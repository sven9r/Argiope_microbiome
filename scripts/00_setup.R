suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(vegan)
  library(broom)
  library(readr)
  library(tibble)
  library(forcats)
  library(ggrepel)
  library(ggpubr)
  library(ggtext)
  library(patchwork)
  library(rstatix)
})

init_dirs <- function() {
  dirs <- c(
    here("data", "processed"),
    here("output", "statistics"),
    here("output", "tables"),
    here("output", "logs"),
    here("figures", "main"),
    here("figures", "supplement"),
    here("figures", "exploratory")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

log_line <- function(x) {
  cat(
    sprintf("- [%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), x),
    file = here("output", "logs", "analysis_log.md"),
    append = TRUE
  )
}

write_session <- function() {
  con <- file(here("output", "logs", "session_info.txt"), "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c("# Session Info", paste("Generated:", Sys.time()), ""), con)
  capture.output(sessionInfo(), file = con, append = TRUE)
}

init_dirs()
log_line("00_setup completed")

allowed_main_figures <- c(
  "Fig01_map.png",
  "fig02_host_phylogeny.png",
  "fig03_historical_dusa_timeline.png",
  "fig04_bacteria_prey_composition.png",
  "fig05_dominant_bacteria_nmds_boxplot.png"
)

allowed_supplement_figures <- c(
  "S01_rarefaction_curves.png",
  "S02_museum_success.png",
  "S03_dna_quality.png",
  "S04_taxon_and_prey_stack.png",
  "S05_logratio_lollipop.png",
  "S06_betadispar.png",
  "S07_symbiont_haplotypes.png",
  "S08_symbiont_16S_distance.png",
  "S09_mantel.png"
)

save_main_figure <- function(plot, filename, width, height, dpi = 300) {
  if (!filename %in% allowed_main_figures) {
    stop("Refusing to write non-whitelisted main figure: ", filename)
  }
  ggplot2::ggsave(
    here::here("figures", "main", filename),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

save_supplement_figure <- function(plot, filename, width, height, dpi = 300) {
  if (!filename %in% allowed_supplement_figures) {
    stop("Refusing to write non-whitelisted supplementary figure: ", filename)
  }
  ggplot2::ggsave(
    here::here("figures", "supplement", filename),
    plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

load_argiope_raw <- function() {
  data_path <- here::here("data", "Raw_data_COI_16S_2.xlsm")
  if (!file.exists(data_path)) stop("Cannot find primary xlsm data file")

  list(
    prey = readxl::read_excel(data_path, sheet = "COI_raw_data"),
    prey_taxa = readxl::read_excel(data_path, sheet = "COI_taxa"),
    bacteria = readxl::read_excel(data_path, sheet = "16S_raw_data"),
    bacteria_taxa = readxl::read_excel(data_path, sheet = "16S_taxa")
  )
}

candidate_terms <- c(
  "dusa",
  "unknown symbiont",
  "mesoplasma",
  "spiroplasma",
  "mycoplasmatota",
  "mollicutes",
  "mycoplasma",
  "tenericutes"
)

tag_candidate_features <- function(bacteria_taxa) {
  tax_long <- bacteria_taxa %>% mutate(across(everything(), as.character))
  tax_text <- apply(tax_long %>% select(-OTU), 1, function(row) paste(row, collapse = " | "))
  hits <- matrix(FALSE, nrow = nrow(tax_long), ncol = length(candidate_terms))
  colnames(hits) <- candidate_terms
  for (i in seq_along(candidate_terms)) {
    hits[, i] <- stringr::str_detect(stringr::str_to_lower(tax_text), candidate_terms[i])
  }
  bind_cols(tax_long, as_tibble(hits), tax_text = tax_text) %>%
    filter(if_any(all_of(candidate_terms), identity))
}

make_dominance_state <- function(bacteria, bacteria_taxa) {
  bacteria_otu <- grep("^OTU", names(bacteria), value = TRUE)
  candidates <- tag_candidate_features(bacteria_taxa)
  otu_dusa <- candidates %>% filter(`dusa` | `unknown symbiont`) %>% pull(OTU) %>% unique()
  otu_meso <- candidates %>% filter(`mesoplasma`) %>% pull(OTU) %>% unique()

  cont_ab <- bacteria %>%
    filter(Status == "Contemporary", str_detect(str_to_lower(Taxon), "bruennichi"))

  cont_long <- cont_ab %>%
    select(ID, all_of(bacteria_otu)) %>%
    pivot_longer(-ID, names_to = "OTU", values_to = "reads")

  state_tbl <- cont_long %>%
    group_by(ID) %>%
    summarize(
      total_16S_reads = sum(reads, na.rm = TRUE),
      dusa_reads = sum(reads[OTU %in% otu_dusa], na.rm = TRUE),
      mesoplasma_reads = sum(reads[OTU %in% otu_meso], na.rm = TRUE),
      other_reads = total_16S_reads - dusa_reads - mesoplasma_reads,
      .groups = "drop"
    ) %>%
    mutate(
      dusa_rel = dusa_reads / pmax(total_16S_reads, 1),
      mesoplasma_rel = mesoplasma_reads / pmax(total_16S_reads, 1),
      other_rel = pmax(other_reads / pmax(total_16S_reads, 1), 0),
      focal_reads = dusa_reads + mesoplasma_reads,
      dusa_focal_rel = if_else(focal_reads > 0, dusa_reads / focal_reads, NA_real_),
      mesoplasma_focal_rel = if_else(focal_reads > 0, mesoplasma_reads / focal_reads, NA_real_),
      log10_dusa_mesoplasma = log10((dusa_reads + 1) / (mesoplasma_reads + 1))
    )

  # Retain the largest overall feature as descriptive metadata only.
  # Focal grouping uses DUSA versus Mesoplasma; other bacteria never set the group.
  largest_feature <- cont_long %>% group_by(ID) %>%
    slice_max(reads, n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(ID, largest_bacterial_feature = OTU, largest_bacterial_feature_reads = reads)
  state_tbl <- state_tbl %>% left_join(largest_feature, by = "ID") %>%
    mutate(
      ambiguous = focal_reads <= 0 |
        (dusa_focal_rel >= 0.2 & mesoplasma_focal_rel >= 0.2),
      max_group = if_else(dusa_reads >= mesoplasma_reads, "DUSA_dominant", "Mesoplasma_dominant"),
      dominance_state = if_else(ambiguous, "ambiguous", max_group),
      profile_group_detail = case_when(
        ambiguous ~ "Mixed", dominance_state == "DUSA_dominant" ~ "DUSA-like",
        TRUE ~ "Mesoplasma"),
      is_TH81_or_TH110 = ID %in% c("TH81", "TH110")
    )

  list(
    state_tbl = state_tbl,
    cont_ab = cont_ab,
    otu_dusa = otu_dusa,
    otu_meso = otu_meso,
    bacteria_otu = bacteria_otu
  )
}

sample_display_label <- function(id) {
  id_chr <- as.character(id)
  th_num <- suppressWarnings(as.integer(stringr::str_match(id_chr, "^TH(\\d+)$")[, 2]))
  dplyr::if_else(
    !is.na(th_num) & th_num >= 81 & th_num <= 110,
    paste0("AG", th_num - 80),
    id_chr
  )
}

# Order specimens by focal DUSA fraction, retaining the previous display identities.
build_focal_arg_map <- function(state_tbl, previous_map) {
  stopifnot(!anyDuplicated(previous_map$ID), !anyDuplicated(state_tbl$ID),
            setequal(previous_map$ID, state_tbl$ID))
  previous <- previous_map %>% transmute(ID, previous_arg_label =
    if ("previous_arg_label" %in% names(previous_map)) previous_arg_label else arg_label)
  state_tbl %>% mutate(focal_order_fraction = dusa_reads / (dusa_reads + mesoplasma_reads)) %>%
    arrange(desc(focal_order_fraction), desc(dusa_reads), ID) %>%
    left_join(previous, by = "ID") %>%
    mutate(arg_label = paste0("Arg", row_number()), arg_order = row_number(),
      ordering_rule = "Descending DUSA/(DUSA+Mesoplasma), then DUSA reads descending, then specimen ID") %>%
    select(ID, arg_label, previous_arg_label, arg_order, ordering_rule)
}

# Arg labels follow focal DUSA fraction and stay fixed across marker filters.
make_fig04_arg_map <- function(state_tbl, cont_ab, min_bacteria_reads = NULL) {
  label_path <- here::here("data", "reference", "arg_sample_map.csv")
  labels <- readr::read_csv(label_path, show_col_types = FALSE)
  stopifnot(nrow(labels) == 30L, !anyDuplicated(labels$ID), !anyDuplicated(labels$arg_label),
            setequal(labels$ID, cont_ab$ID),
            identical(sort(as.integer(sub("^Arg", "", labels$arg_label))), 1:30))
  expected <- build_focal_arg_map(state_tbl, labels)
  stopifnot(identical(labels$ID[order(as.integer(sub("^Arg", "", labels$arg_label)))], expected$ID))
  result <- state_tbl %>%
    inner_join(cont_ab %>% select(ID, Number, BacteriaReads, DominantSymbiont, UnknownSymbiont), by = "ID") %>%
    inner_join(labels, by = "ID") %>% arrange(as.integer(sub("^Arg", "", arg_label)))
  if (!is.null(min_bacteria_reads)) result <- result %>% filter(BacteriaReads >= min_bacteria_reads)
  result
}

contemporary_analysis_cohort <- function(state, prey) {
  result <- make_fig04_arg_map(state$state_tbl, state$cont_ab) %>%
    left_join(prey %>% select(ID, PreyReads), by = "ID") %>%
    mutate(eligible_16S = BacteriaReads >= 1500, eligible_COI = PreyReads >= 1693,
           paired = eligible_16S & eligible_COI,
           clear_paired = paired & dominance_state != "ambiguous")
  stopifnot(nrow(result) == 30L, sum(result$eligible_16S) == 29L,
            sum(result$eligible_COI) == 29L, sum(result$paired) == 28L,
            sum(result$clear_paired) == 26L,
            setequal(result$ID[!result$paired], c("TH84", "TH99")),
            setequal(result$ID[result$paired & !result$clear_paired], c("TH81", "TH110")),
            sum(result$paired & result$dominance_state == "DUSA_dominant") == 16L,
            sum(result$paired & result$dominance_state == "Mesoplasma_dominant") == 10L)
  result
}

apply_fig04_arg_labels <- function(data, fig04_arg_map, fallback = c("next_arg", "ID", "sample_display_label")) {
  fallback <- match.arg(fallback)
  next_arg_start <- max(as.integer(stringr::str_remove(fig04_arg_map$arg_label, "^Arg")), na.rm = TRUE)

  data %>%
    left_join(fig04_arg_map %>% select(ID, fig04_arg_label = arg_label), by = "ID") %>%
    mutate(
      arg_label = case_when(
        !is.na(fig04_arg_label) ~ fig04_arg_label,
        fallback == "next_arg" ~ paste0("Arg", next_arg_start + cumsum(is.na(fig04_arg_label))),
        fallback == "sample_display_label" ~ sample_display_label(ID),
        TRUE ~ as.character(ID)
      )
    ) %>%
    select(-fig04_arg_label)
}

write_state_tables <- function(state) {
  readr::write_csv(state$state_tbl, here::here("output", "tables", "table_3_dominance_state_per_sample.csv"))
  readr::write_csv(
    state$state_tbl %>% select(ID, log10_dusa_mesoplasma),
    here::here("output", "tables", "table_3_dusa_mesoplasma_logratio.csv")
  )
  readr::write_csv(
    state$state_tbl %>% count(dominance_state, name = "n"),
    here::here("output", "tables", "table_3_dominance_state_counts.csv")
  )
  rho <- suppressWarnings(cor(
    state$state_tbl$dusa_rel,
    state$state_tbl$mesoplasma_rel,
    method = "spearman",
    use = "complete.obs"
  ))
  readr::write_csv(
    tibble(statistic = "spearman_rho_dusa_vs_mesoplasma", value = rho),
    here::here("output", "statistics", "03_dusa_mesoplasma_exclusion.csv")
  )
}

prey_order_palette <- c(
  "Blattodea" = "#6a3d9a",
  "Coleoptera" = "#a50013",
  "Diptera" = "#fb6449",
  "Hemiptera" = "#eeee00",
  "Hymenoptera" = "#74c476",
  "Lepidoptera" = "#6baed6",
  "Orthoptera" = "#2878b8",
  "Other (<3%)" = "#999999",
  "Other" = "#999999",
  "Unknown" = "#999999"
)

symbiont_component_palette <- c(
  "Other (<10%)" = "#999999",
  "*Mesoplasma sp.*" = "#eeee00",
  "DUSA" = "#ff6a0a"
)

dominance_palette <- c(
  "Ambiguous" = "black",
  "DUSA" = "#ff6a0a",
  "*Mesoplasma sp.*" = "#eeee00"
)

log_line("00_setup completed")
