#!/usr/bin/env Rscript
# Supplementary Figure 3: read recovery and thresholded OTU counts by collection period.
source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
periods <- c("veryold", "old", "medium", "young", "contemporary")
period_labels <- c("\u22641900", "1901\u20131950", "1951\u20132000", "2001\u20132020", "Contemporary")
original_letters <- list(A = c("b", "b", "b", "b", "a"),
                         B = c("b", "b", "b", "a", "a"),
                         C = c("ab", "a", "ab", "b", "a"))
expected_n <- list(A = c(10, 20, 29, 18, 30), B = c(10, 20, 29, 18, 30),
                   C = c(1, 3, 12, 10, 29), D = c(2, 2, 10, 15, 29))
panel_data <- list()

for (marker in c("16S", "COI")) {
  data <- if (marker == "16S") raw$bacteria else raw$prey
  depth <- if (marker == "16S") "BacteriaReads" else "PreyReads"
  profile_threshold <- if (marker == "16S") 507 else 1133
  feature_threshold <- if (marker == "16S") 20 else 100
  otus <- grep("^OTU", names(data), value = TRUE)
  # TotalOTU is the workbook's count of features meeting this per-specimen
  # read threshold. It is not richness recalculated after rarefaction.
  stopifnot(nrow(data) == 112L, !anyDuplicated(data$ID),
            all(data$TotalOTU == rowSums(data[, otus] >= feature_threshold)))
  data <- data %>%
    mutate(period = cut(suppressWarnings(as.numeric(Year)),
                        breaks = c(-Inf, 1900, 1950, 2000, 2020, Inf), labels = periods)) %>%
    filter(!is.na(period))
  for (metric in c(depth, "TotalOTU")) {
    panel <- if (marker == "16S") {
      if (metric == depth) "A" else "C"
    } else {
      if (metric == depth) "B" else "D"
    }
    selected <- if (metric == "TotalOTU") filter(data, .data[[depth]] >= profile_threshold) else data
    panel_data[[panel]] <- selected %>%
      transmute(panel, marker, ID, Status, Year, period, metric, value = .data[[metric]])
    stopifnot(identical(as.integer(table(selected$period)), as.integer(expected_n[[panel]])))
  }
}
panel_data <- bind_rows(panel_data) %>% mutate(panel = factor(panel, levels = c("A", "B", "C", "D")))
write_csv(panel_data, here::here("output", "tables", "S03_panel_data.csv"))

kruskal <- list()
dunn <- list()
panels <- list()
for (panel in c("A", "B", "C", "D")) {
  data <- filter(panel_data, .data$panel == .env$panel)
  counts <- data %>% count(period, name = "n") %>%
    mutate(label = paste0(period_labels, "<br>n=", n))
  test <- kruskal.test(value ~ period, data = data)
  kruskal[[panel]] <- tibble(panel, marker = first(data$marker), metric = first(data$metric),
                            n = nrow(data), chi_squared = unname(test$statistic),
                            df = unname(test$parameter), P = test$p.value)
  data <- left_join(data, counts, by = "period") %>%
    mutate(label = factor(label, levels = counts$label))
  annotations <- data %>% group_by(period, label) %>% summarize(max = max(value), .groups = "drop")
  if (panel %in% names(original_letters)) {
    # Dunn's pooled-rank, tie-corrected two-sided comparisons with Bonferroni
    # adjustment reproduce the original FSA dunnTest(method = "bonferroni").
    comparisons <- rstatix::dunn_test(data, value ~ period, p.adjust.method = "bonferroni")
    letters <- setNames(original_letters[[panel]], periods)
    shared <- mapply(function(a, b) {
      length(intersect(strsplit(letters[[a]], "")[[1]], strsplit(letters[[b]], "")[[1]])) > 0
    }, comparisons$group1, comparisons$group2)
    # Keep the original letter convention only while every comparison agrees.
    stopifnot(all(shared == (comparisons$p.adj >= 0.05)))
    # Report Z in group1-minus-group2 order.
    dunn[[panel]] <- comparisons %>%
      transmute(panel, group1, group2, n1, n2, Z = -statistic,
                P_unadjusted = p, P_bonferroni = p.adj,
                letter1 = letters[group1], letter2 = letters[group2])
    annotations$letter <- letters[as.character(annotations$period)]
  }
  set.seed(100)
  p <- ggplot(data, aes(label, value, fill = period)) +
    geom_jitter(show.legend = FALSE,
                position = position_jitterdodge(jitter.width = if (panel == "D") 0.4 else 1.5,
                                                dodge.width = 0.7), shape = 21, size = 3.8) +
    stat_summary(fun = median, geom = "crossbar", position = position_dodge(width = 0.8),
                 width = 0.6, show.legend = FALSE, linewidth = 0.35) +
    ggpubr::stat_compare_means(method = "kruskal.test") +
    scale_fill_brewer(palette = if (panel %in% c("A", "C")) "Blues" else "Greens") +
    scale_y_continuous(limits = c(0, if (panel == "A") 10000 else if (panel == "B") 20000 else 25)) +
    labs(x = if (panel %in% c("C", "D")) "Year of collection" else NULL,
         y = if (panel == "A") "Reads [n]" else if (panel == "C") "Thresholded OTUs [n]" else NULL,
         title = if (panel == "A") "Bacterial DNA\n" else if (panel == "B") "Prey DNA\n" else NULL) +
    theme_classic() +
    theme(plot.title = element_text(size = 13), legend.position = "none",
          axis.text.x = ggtext::element_markdown(size = 9), axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12), plot.tag = element_text(size = 13, face = "bold"))
  if (panel %in% names(original_letters)) {
    p <- p + geom_text(data = annotations, aes(x = label, y = max, label = letter),
                       inherit.aes = FALSE, vjust = -3)
  }
  panels[[panel]] <- p
}
write_csv(bind_rows(kruskal), here::here("output", "statistics", "S03_kruskal.csv"))
write_csv(bind_rows(dunn), here::here("output", "statistics", "S03_dunn_bonferroni.csv"))
figure <- ((panels$A | panels$B) / (panels$C | panels$D)) +
  plot_annotation(tag_levels = "A")
save_supplement_figure(figure, "S03_dna_quality.png", width = 10.5, height = 8)
