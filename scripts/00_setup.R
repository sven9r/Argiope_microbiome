suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(vegan)
  library(readr)
  library(tibble)
  library(forcats)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
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
      log10_dusa_mesoplasma = log10((dusa_reads + 1) / (mesoplasma_reads + 1))
    )

  max_feat <- cont_long %>%
    group_by(ID) %>%
    slice_max(reads, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      max_group = case_when(
        OTU %in% otu_dusa ~ "DUSA_dominant",
        OTU %in% otu_meso ~ "Mesoplasma_dominant",
        TRUE ~ "Other_dominant"
      )
    ) %>%
    select(ID, max_group)

  state_tbl <- state_tbl %>%
    left_join(max_feat, by = "ID") %>%
    mutate(
      ambiguous = (max_group == "Other_dominant") |
        (dusa_rel >= 0.2 & mesoplasma_rel >= 0.2) |
        abs(log10_dusa_mesoplasma) < 0.5,
      dominance_state = case_when(
        ambiguous ~ "ambiguous",
        max_group == "DUSA_dominant" ~ "DUSA_dominant",
        max_group == "Mesoplasma_dominant" ~ "Mesoplasma_dominant",
        TRUE ~ "ambiguous"
      ),
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

# The deposited mapping keeps Arg labels stable when read thresholds change.
make_fig04_arg_map <- function(state_tbl, cont_ab, min_bacteria_reads = 1500) {
  labels <- readr::read_csv(here::here("data", "reference", "arg_sample_map.csv"), show_col_types = FALSE) %>%
    transmute(ID, arg_label = union_arg_label)
  state_tbl %>%
    inner_join(cont_ab %>% select(ID, Number, BacteriaReads, DominantSymbiont, UnknownSymbiont), by = "ID") %>%
    inner_join(labels, by = "ID") %>%
    filter(BacteriaReads >= min_bacteria_reads) %>%
    arrange(as.integer(sub("^Arg", "", arg_label)))
}

apply_fig04_arg_labels <- function(data, fig04_arg_map, fallback = "next_arg") {
  labels <- readr::read_csv(here::here("data", "reference", "arg_sample_map.csv"), show_col_types = FALSE) %>%
    transmute(ID, arg_label = union_arg_label)
  result <- left_join(data, labels, by = "ID")
  stopifnot(!anyNA(result$arg_label))
  result
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
