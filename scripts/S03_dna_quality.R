#!/usr/bin/env Rscript
# Maintained S3 renderer, preserving the data flow and styling of
# scripts/_archive_20260601/script_analysis_thesis.Rmd:375-683.
# Approved changes: shared A/B read axis; C/D descriptive, without P values/letters.
# Run from the project root with the installed R 4.5 project environment.
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(ggpubr)
  library(patchwork)
})

root <- normalizePath(getwd())
workbook <- file.path(root, "data/Raw_data_COI_16S_2.xlsm")
target <- file.path(root, "figures/supplement/S03_dna_quality.png")
support <- file.path(root, "output/tables")
statistics <- file.path(root, "output/statistics")
dir.create(support, recursive = TRUE, showWarnings = FALSE)
dir.create(statistics, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
prior_sha256 <- if (file.exists(target)) digest::digest(file = target, algo = "sha256") else NA_character_
bacteria <- read_excel(workbook, sheet = "16S_raw_data")
prey <- read_excel(workbook, sheet = "COI_raw_data")
stopifnot(nrow(bacteria) == 112L, nrow(prey) == 112L,
          !anyDuplicated(bacteria$ID), !anyDuplicated(prey$ID))

periods <- c("veryold", "old", "medium", "young", "contemporary")
period_labels <- c("<1900", "1900-1950", "1951-2000", "2001-2020", "Contemporary")
period_cut <- function(y, labels = period_labels)
  cut(y, breaks = c(-Inf, 1900, 1950, 2000, 2020, Inf), labels = labels)
prepare_panel <- function(x, metric) {
  x <- x %>% select(ID, Taxon, Status, all_of(metric), Year) %>% filter(Year != "NA")
  x$Year <- as.numeric(x$Year)
  counts <- x %>% mutate(category = period_cut(Year)) %>% count(category) %>%
    unite("label", category, n, sep = "<br>n=", remove = FALSE)
  x %>% mutate(category = period_cut(Year), Group = period_cut(Year, periods)) %>%
    left_join(counts, by = "category")
}

# Keep the original actual rarefaction and row.names merge, including its order.
# The plotted TotalOTU is retained metadata, not a recount of rarefied OTUs.
rarefied_metadata <- function(x, read_column, depth) {
  set.seed(555)
  m <- GUniFrac::Rarefy(x %>% select(starts_with("OTU")), depth)$otu.tab.rff %>% as.data.frame()
  z <- x %>% select(-starts_with("OTU")) %>% filter(.data[[read_column]] >= depth) %>%
    merge(m, by = "row.names") %>% as.data.frame()
  otus <- grep("^OTU", names(z), value = TRUE)
  stopifnot(all(rowSums(z[, otus, drop = FALSE]) == depth),
    identical(as.numeric(z$TotalOTU), as.numeric(x$TotalOTU[match(z$ID, x$ID)])))
  z
}
panels <- list(A = prepare_panel(bacteria, "BacteriaReads"),
               B = prepare_panel(prey, "PreyReads"),
               C = prepare_panel(rarefied_metadata(bacteria, "BacteriaReads", 507), "TotalOTU"),
               D = prepare_panel(rarefied_metadata(prey, "PreyReads", 1133), "TotalOTU"))
metrics <- c(A = "BacteriaReads", B = "PreyReads", C = "TotalOTU", D = "TotalOTU")
expected_n <- list(A = c(10,20,29,18,30), B = c(10,20,29,18,30),
                   C = c(1,3,12,10,29), D = c(2,2,10,15,29))
expected_medians <- list(A = c(219.5,134.5,342,716.5,5004.5),
  B = c(48,3.5,4,4481,7239.5), C = c(5,2,5.5,9.5,3), D = c(4.5,3,6,8,6))
for (p in names(panels)) {
  z <- panels[[p]]
  stopifnot(identical(as.integer(table(z$Group)), as.integer(expected_n[[p]])),
    identical(as.numeric(tapply(z[[metrics[p]]], z$Group, median)), expected_medians[[p]]))
}

# Preserve the original A/B letters and verify their complete pairwise meaning.
# rstatix is installed; FSA/rcompanion from the thesis runtime are unavailable.
letters_by_panel <- list(A = setNames(c("b","b","b","b","a"), periods),
                        B = setNames(c("b","b","b","a","a"), periods))
original_dunn <- list(
  A = c(1,1,0.727616911773515,9.96184774141458e-7,0.584890860525442,
        0.152583248431988,2.25413523030728e-11,1,1.4481385972752e-7,0.000330935034299017),
  B = c(1,1,0.0194985084731574,0.000135677411405074,1,
        0.000127997197923042,6.33321605052317e-9,0.00970007830297885,1.88592465245423e-6,1))
original_kw <- c(A = 3.24135260326883e-13, B = 2.16073218572445e-11)
original_kw_H <- c(A = 64.5240641286608, B = 55.8452252328997)
kw_proof <- list(); dunn_proof <- list(); annotation <- list()
pair_names <- combn(periods, 2, simplify = FALSE)
for (p in c("A", "B")) {
  z <- panels[[p]]; metric <- metrics[p]; f <- reformulate("Group", response = metric)
  kw <- kruskal.test(f, data = z)
  stopifnot(isTRUE(all.equal(unname(kw$p.value), unname(original_kw[p]), tolerance = 1e-10)),
    isTRUE(all.equal(unname(kw$statistic), unname(original_kw_H[p]), tolerance = 1e-12)),
    unname(kw$parameter) == 4L)
  d <- rstatix::dunn_test(z, f, p.adjust.method = "bonferroni", detailed = TRUE)
  wanted <- vapply(pair_names, paste, collapse = "/", FUN.VALUE = character(1))
  i <- match(wanted, paste(d$group1, d$group2, sep = "/")); stopifnot(!anyNA(i))
  d <- d[i, ]; d$original_P_adjusted <- original_dunn[[p]]
  d$absolute_P_difference <- abs(d$p.adj - d$original_P_adjusted)
  shared <- vapply(seq_len(nrow(d)), function(i) {
    a <- strsplit(letters_by_panel[[p]][as.character(d$group1[i])], "")[[1]]
    b <- strsplit(letters_by_panel[[p]][as.character(d$group2[i])], "")[[1]]
    length(intersect(a, b)) > 0
  }, logical(1))
  stopifnot(nrow(d) == 10L, all(d$absolute_P_difference < 1e-12),
            all(shared == (d$p.adj >= 0.05)))
  d$panel <- p; d$original_letter_relations_preserved <- shared == (d$p.adj >= 0.05)
  dunn_proof[[p]] <- d %>% select(any_of(c("panel", "group1", "group2", "n1", "n2",
    "statistic", "p", "p.adj", "original_P_adjusted", "absolute_P_difference",
    "original_letter_relations_preserved")))
  kw_proof[[p]] <- data.frame(panel = p, n = nrow(z), H = unname(kw$statistic),
    df = unname(kw$parameter), P = kw$p.value, original_P = original_kw[p])
  annotation[[p]] <- z %>% group_by(Group, category, label) %>%
    summarize(max = max(.data[[metric]]), .groups = "drop") %>%
    mutate(Letter = unname(letters_by_panel[[p]][as.character(Group)]))
}

make_panel <- function(p, read_limits = c(0, 20000)) {
  z <- panels[[p]]; metric <- metrics[p]; reads <- p %in% c("A", "B")
  bacterial <- p %in% c("A", "C")
  # Exactly retain source seeding before the jitter-position constructor. The
  # constructor samples its own seed; passing seed=100 instead changes offsets.
  set.seed(100)
  source_width <- if (p == "D") 0.4 else 1.5
  pos <- position_jitterdodge(jitter.width = source_width, dodge.width = 0.7)
  stopifnot(pos$seed == 1321812471L)
  # ggplot2 3.4.2 divided jitter width by all fill factor levels + 2. A-C use
  # five category levels; D uses the character label and thus has zero levels.
  # Current ggplot2 counts the one group at each x instead.
  # Preserve the original effective width without changing its RNG or row order.
  # https://github.com/tidyverse/ggplot2/blob/v3.4.2/R/position-jitterdodge.R#L45-L56
  probe <- data.frame(x = 1:5, y = 1:5, group = 1:5, PANEL = 1,
                      fill = if (p == "D") periods else factor(periods, levels = periods))
  actual_width <- pos$setup_params(probe)$jitter.width
  legacy_divisor <- if (p == "D") 2 else 7
  pos$jitter.width <- source_width * (source_width / legacy_divisor) / actual_width
  g <- ggplot(z, aes(x = label, y = .data[[metric]],
                     fill = .data[[if (p == "D") "label" else "category"]])) +
    geom_jitter(show.legend = FALSE, position = pos, pch = 21, size = 3.8) +
    stat_summary(fun = median, geom = "crossbar", position = position_dodge(width = 0.8),
      width = 0.6, show.legend = FALSE, linewidth = 0.35) +
    labs(x = if (reads) "" else "Year of catch\n",
      y = if (p == "A") "Reads [n]" else if (p == "C") "zOTUs [n]" else "",
      title = if (p == "A") "Bacterial DNA\n" else if (p == "B") "Prey DNA\n" else NULL) +
    theme_classic() + theme(plot.title = element_text(size = 13),
      legend.title = element_blank(), axis.text.x = element_markdown(size = 11),
      axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 12, angle = 90),
      axis.title.x = element_text(size = 12, vjust = -2)) +
    scale_fill_brewer(palette = if (bacterial) "Blues" else "Greens")
  if (reads) g <- g +
    stat_compare_means(method = "kruskal.test") +
    scale_y_continuous(limits = read_limits, breaks = seq(0, read_limits[2], by = 5000)) +
    geom_text(data = annotation[[p]], aes(x = label, y = max, label = Letter), vjust = -3)
  else g <- g + scale_y_continuous(limits = c(0, 25))
  g
}
plots <- lapply(names(panels), make_panel); names(plots) <- names(panels)
built <- lapply(plots, ggplot_build)
for (p in names(plots)) {
  x <- as.numeric(built[[p]]$data[[1]]$x)
  stopifnot(max(abs(x - round(x))) <= (if (p == "D") 0.4/2 else 1.5/7) + 1e-12)
}
for (p in c("A", "B")) {
  s <- plots[[p]]$scales$get_scales("y")
  stopifnot(identical(s$limits, c(0, 20000)), identical(s$breaks, seq(0,20000,5000)))
  original_scale <- ggplot_build(make_panel(p, if (p == "A") c(0,10000) else c(0,20000)))
  stopifnot(identical(built[[p]]$data[[1]]$x, original_scale$data[[1]]$x),
            identical(built[[p]]$data[[1]]$y, original_scale$data[[1]]$y))
}
stopifnot(!any(vapply(plots$C$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))),
          !any(vapply(plots$D$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))),
          length(plots$C$layers) == 2L, length(plots$D$layers) == 2L)

plot_data <- bind_rows(lapply(names(panels), function(p) {
  z <- panels[[p]]
  data.frame(panel = p, source_row_order = seq_len(nrow(z)), ID = z$ID,
    Status = z$Status, Year = z$Year, period = as.character(z$Group),
    metric = unname(metrics[p]), value = z[[metrics[p]]])
}))
summary <- plot_data %>% group_by(panel, period) %>%
  summarize(n = n(), median = median(value), min = min(value), max = max(value), .groups = "drop")
figure <- ((plots$A | plots$B) / (plots$C | plots$D)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) &
  theme(plot.tag = element_text(size = 13, face = "bold"))
ggsave(target, figure, width = 12, height = 9, dpi = 300, bg = "white")
readr::write_csv(plot_data, file.path(support, "S03_render_plot_data.csv"))
readr::write_csv(summary, file.path(support, "S03_render_group_summary.csv"))
readr::write_csv(bind_rows(kw_proof), file.path(support, "S03_render_AB_kruskal_proof.csv"))
readr::write_csv(bind_rows(dunn_proof), file.path(support, "S03_render_AB_dunn_proof.csv"))
readr::write_csv(plot_data %>% mutate(marker = if_else(panel %in% c("A", "C"), "16S", "COI")),
  file.path(support, "S03_panel_data.csv"))
readr::write_csv(bind_rows(kw_proof) %>% transmute(panel,
  marker = if_else(panel == "A", "16S", "COI"),
  metric = if_else(panel == "A", "BacteriaReads", "PreyReads"), n, chi_squared = H, df, P),
  file.path(statistics, "S03_kruskal.csv"))
readr::write_csv(bind_rows(dunn_proof) %>% rowwise() %>% transmute(panel, group1, group2, n1, n2,
  Z = -statistic, P_unadjusted = p, P_bonferroni = p.adj,
  letter1 = letters_by_panel[[panel]][[as.character(group1)]],
  letter2 = letters_by_panel[[panel]][[as.character(group2)]]) %>% ungroup(),
  file.path(statistics, "S03_dunn_bonferroni.csv"))
audit <- list(figure = "figures/supplement/S03_dna_quality.png",
  source = "scripts/_archive_20260601/script_analysis_thesis.Rmd:375-683",
  approved_changes = c("A/B shared 0-20000 read axis, ticks every 5000",
    "C/D descriptive: no omnibus P values or pairwise letters"),
  width_pixels = 3600L, height_pixels = 2700L, panel_n = vapply(panels, nrow, integer(1)),
  contemporary_n = vapply(panels, function(z) sum(z$Status == "Contemporary"), integer(1)),
  period_Ns_and_medians_preserved = TRUE, AB_tests_and_letters_preserved = TRUE,
  AB_point_coordinates_preserved_across_scale_change = TRUE,
  jitter_constructor_seed = 1321812471L,
  jitter_effective_width = c(A = 1.5/7, B = 1.5/7, C = 1.5/7, D = 0.4/2),
  jitter_compatibility = "Retain ggplot2 3.4.2 width/(fill factor levels+2): five levels in A-C, character fill in D; compensate for changed current dodge counting without changing RNG or row order",
  CD_plotted_quantity = "Imported TotalOTU metadata retained in the original rarefied-profile data flow",
  prior_output_sha256 = prior_sha256,
  workbook_sha256 = digest::digest(file = workbook, algo = "sha256"),
  output_sha256 = digest::digest(file = target, algo = "sha256"))
jsonlite::write_json(audit, file.path(support, "S03_render_audit.json"), pretty = TRUE, auto_unbox = TRUE)
cat("S3 updated:", target, "\nPanel contemporary N:",
    paste(audit$contemporary_n, collapse = "/"), "\n")
