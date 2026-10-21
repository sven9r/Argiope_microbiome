# Canonical Figure 5: 1,500/1,693 paired cohort and focal-only DUSA/Mesoplasma grouping,
# prey rarefaction to 1,693 reads, BH adjustment and deterministic seed 100.
# Produces full model evidence under output/figure5 and publication artwork.
root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "output", "figure5")
opt <- list(`state-rule` = "focal-only", `order-transform` = "rarefied", `p-adjust` = "BH", seed = "100")
seed <- 100L
source(file.path(root, "scripts/00_setup.R"))
workbook <- file.path(root, "data/Raw_data_COI_16S_2.xlsm")
map_path <- file.path(root, "data/reference/arg_sample_map.csv")
stopifnot(file.exists(workbook), file.exists(map_path))
bac <- read_excel(workbook, sheet = "16S_raw_data")
prey <- read_excel(workbook, sheet = "COI_raw_data")
taxa <- read_excel(workbook, sheet = "COI_taxa")
bacterial_taxa <- read_excel(workbook, sheet = "16S_taxa")
is_focal <- function(d) d$Status == "Contemporary" & grepl("bruennichi", tolower(d$Taxon))
bac <- bac[is_focal(bac), , drop = FALSE]
prey <- prey[is_focal(prey), , drop = FALSE]
stopifnot(nrow(bac) == 30, nrow(prey) == 30, !anyDuplicated(bac$ID), !anyDuplicated(prey$ID), setequal(bac$ID, prey$ID))
bot <- grep("^OTU", names(bac), value = TRUE)
pot <- grep("^OTU", names(prey), value = TRUE)
bm <- as.matrix(bac[, bot]); rownames(bm) <- bac$ID
pm <- as.matrix(prey[, pot]); rownames(pm) <- prey$ID
stopifnot(all(is.finite(bm)), all(is.finite(pm)), all(bm >= 0), all(pm >= 0),
          all(bm == floor(bm)), all(pm == floor(pm)),
          all(rowSums(bm) == bac$BacteriaReads), all(rowSums(pm) == prey$PreyReads),
          all(c("OTU2", "OTU3", "OTU145") %in% bot))
# Use the shared focal-only rule; other bacterial reads never determine this group.
state <- make_dominance_state(bac, bacterial_taxa)$state_tbl %>%
  transmute(ID, bacterial_reads = total_16S_reads, dusa_reads, mesoplasma_reads,
    focal_reads = dusa_reads + mesoplasma_reads,
    dusa_rel, mesoplasma_rel,
    dusa_fraction_focal = dusa_reads / focal_reads,
    mesoplasma_fraction_focal = mesoplasma_reads / focal_reads,
    log10_ratio = log10_dusa_mesoplasma, ambiguous,
    strict_group = recode(dominance_state, DUSA_dominant = "DUSA-like",
      Mesoplasma_dominant = "Mesoplasma", ambiguous = "Mixed"))
labels <- read_csv(map_path, show_col_types = FALSE) %>% select(ID, figure4_label = arg_label)
stopifnot(!anyDuplicated(labels$ID), !anyDuplicated(labels$figure4_label), setequal(labels$ID, bac$ID))
cohort <- state %>% left_join(prey %>% select(ID, prey_reads = PreyReads), by = "ID") %>%
  left_join(labels, by = "ID") %>% mutate(
    eligible_16S = bacterial_reads >= 1500, eligible_COI = prey_reads >= 1693,
    paired = eligible_16S & eligible_COI,
    display_group = strict_group,
    in_model = paired & !ambiguous,
    model_group = ifelse(in_model, display_group, NA_character_)
  ) %>% arrange(as.integer(sub("Arg", "", figure4_label)))
stopifnot(sum(cohort$paired) == 28, setequal(cohort$ID[!cohort$paired], c("TH84", "TH99")),
          all(cohort$strict_group %in% c("DUSA-like", "Mesoplasma", "Mixed")))
stopifnot(sum(cohort$in_model) == 26, setequal(cohort$ID[cohort$paired & cohort$ambiguous], c("TH81", "TH110")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv(cohort, file.path(out, "cohort_all30.csv"))
write_csv(cohort %>% filter(paired), file.path(out, "display_cohort.csv"))
write_csv(cohort %>% filter(in_model), file.path(out, "model_cohort.csv"))
source_files <- c(workbook = workbook, figure4_map = map_path,
                  reference_state_helper = file.path(root, "scripts/00_setup.R"),
                  canonical_source = file.path(root, "scripts/05_community_comparison.R"))
stopifnot(all(file.exists(source_files)))
write_csv(tibble(source = names(source_files), path = source_files, md5 = unname(tools::md5sum(source_files))), file.path(out, "input_hashes.csv"))
suppressPackageStartupMessages({library(vegan); library(ggplot2); library(patchwork); library(ggrepel)})
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seeds <- c(distance_16S = seed, distance_COI = seed + 1L, nmds_16S = seed + 2L,
           nmds_COI = seed + 3L, permanova_16S = seed + 4L, permanova_COI = seed + 5L,
           dispersion_16S = seed + 6L, dispersion_COI = seed + 7L, order_rarefaction = seed + 8L)
display <- cohort %>% filter(paired) %>% arrange(ID)
# Statistical input order follows immutable specimen IDs, independently of Arg display labels.
ids <- display$ID
model <- cohort %>% filter(in_model) %>% arrange(ID) %>% mutate(model_group = factor(model_group, levels = c("DUSA-like", "Mesoplasma")))
pal <- c(`DUSA-like` = "#ff6a0a", Mesoplasma = "#eeee00", Mixed = "black")
analyse_marker <- function(mat, depth, marker, dist_seed, ord_seed, perm_seed, bd_seed) {
  mat <- mat[ids, , drop = FALSE]
  stopifnot(all(rowSums(mat) >= depth), identical(rownames(mat), display$ID))
  set.seed(dist_seed)
  full_distance <- avgdist(mat, sample = depth, dmethod = "bray", iterations = 100)
  stopifnot(attr(full_distance, "Size") == 28L, identical(attr(full_distance, "Labels"), ids))
  set.seed(ord_seed)
  ord <- metaMDS(full_distance, k = 2, trymax = 100, trace = FALSE)
  scores <- as.data.frame(vegan::scores(ord, display = "sites")) %>%
    tibble::rownames_to_column("ID") %>% left_join(display, by = "ID")
  test_distance <- as.dist(as.matrix(full_distance)[model$ID, model$ID, drop = FALSE])
  stopifnot(attr(test_distance, "Size") == nrow(model), identical(attr(test_distance, "Labels"), model$ID))
  set.seed(perm_seed)
  ad <- adonis2(test_distance ~ model_group, data = model, permutations = 999)
  bd <- betadisper(test_distance, model$model_group, type = "median", bias.adjust = FALSE)
  set.seed(bd_seed)
  bd_test <- permutest(bd, permutations = 999)
  stats <- tibble(marker = marker, n_display = length(ids), n_model = attr(test_distance, "Size"),
                  group1 = "DUSA-like", n_group1 = sum(model$model_group == "DUSA-like"),
                  group2 = "Mesoplasma", n_group2 = sum(model$model_group == "Mesoplasma"),
                  F = ad$F[1], R2 = ad$R2[1], p = ad$`Pr(>F)`[1], stress = ord$stress,
                  dispersion_F = bd_test$tab[1, "F"], dispersion_p = bd_test$tab[1, "Pr(>F)"],
                  dispersion_center = "spatial median", bias_adjust = FALSE)
  write_csv(scores, file.path(out, paste0(marker, "_ordination_scores.csv")))
  write_csv(stats, file.path(out, paste0(marker, "_statistics.csv")))
  write_csv(tibble(ID = model$ID, group = model$model_group, distance_to_spatial_median = bd$distances), file.path(out, paste0(marker, "_dispersion_values.csv")))
  saveRDS(list(matrix = mat, full_distance = full_distance, model_distance = test_distance,
               model_metadata = model, ordination = ord, permanova = ad, dispersion = bd,
               dispersion_test = bd_test, statistics = stats), file.path(out, paste0(marker, "_objects.rds")))
  shown_p <- format(stats$p, scientific = FALSE, trim = TRUE)
  displayed_counts <- table(factor(display$display_group, levels = c("DUSA-like", "Mesoplasma", "Mixed")))
  displayed_labels <- c(
    `DUSA-like` = paste0("DUSA-like (n=", displayed_counts[["DUSA-like"]], ")"),
    Mesoplasma = paste0("Mesoplasma (n=", displayed_counts[["Mesoplasma"]], ")"),
    Mixed = paste0("Mixed (n=", displayed_counts[["Mixed"]], ")")
  )
  p <- ggplot(scores, aes(NMDS1, NMDS2, fill = display_group)) +
    geom_point(shape = 21, color = "black", size = 3.5, stroke = .4) +
    geom_text_repel(aes(label = ifelse(ambiguous, figure4_label, "")), seed = seed + 9L, size = 3) +
    scale_fill_manual(values = pal, drop = FALSE, labels = displayed_labels) + theme_classic(base_size = 12) +
    labs(x = sprintf("MDS1\nStress = %.1f%%", 100 * ord$stress), y = "MDS2", fill = "Focal read profile",
         title = ifelse(marker == "16S", "Bacterial read profiles", "Prey read profiles"))
  list(plot = p, statistics = stats)
}
a <- analyse_marker(bm, 1500, "16S", seeds[1], seeds[3], seeds[5], seeds[7])
b <- analyse_marker(pm, 1693, "COI", seeds[2], seeds[4], seeds[6], seeds[8])

# Panel C uses the same tested cohort, with zeros retained for every order.
order_mat <- pm[model$ID, , drop = FALSE]
set.seed(seeds[["order_rarefaction"]])
order_mat <- vegan::rrarefy(order_mat, sample = 1693)
stopifnot(all(rowSums(order_mat) == 1693))
stopifnot(!anyDuplicated(taxa$OTU), all(pot %in% taxa$OTU))
major <- c("Blattodea", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Orthoptera")
order_names <- c(major, "Other prey orders")
order_map <- taxa$Order[match(pot, taxa$OTU)]
order_map[is.na(order_map) | !order_map %in% major] <- "Other prey orders"
aggregated <- sapply(order_names, function(x) rowSums(order_mat[, order_map == x, drop = FALSE]))
rownames(aggregated) <- model$ID
proportions <- 100 * aggregated / rowSums(aggregated)
stopifnot(all(abs(rowSums(proportions) - 100) < 1e-8))
order_data <- as.data.frame(proportions, check.names = FALSE) %>% tibble::rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "Order", values_to = "percent") %>%
  left_join(model %>% select(ID, figure4_label, model_group), by = "ID") %>%
  mutate(Order = factor(Order, levels = order_names))
tests <- bind_rows(lapply(order_names, function(o) {
  d <- order_data %>% filter(Order == o)
  t <- wilcox.test(d$percent[d$model_group == "DUSA-like"],
                   d$percent[d$model_group == "Mesoplasma"],
                   paired = FALSE, exact = FALSE, correct = TRUE)
  tibble(Order = o, W = unname(t$statistic), p_raw = t$p.value,
         median_difference_DUSA_minus_Meso = median(d$percent[d$model_group == "DUSA-like"]) - median(d$percent[d$model_group == "Mesoplasma"]))
})) %>% mutate(p_adjusted = p.adjust(p_raw, method = opt$`p-adjust`), adjustment = opt$`p-adjust`,
              label = paste0("BH P = ", signif(p_adjusted, 3)), x = seq_along(order_names))
write_csv(order_data, file.path(out, "panelC_values.csv"))
write_csv(tests, file.path(out, "panelC_WMW_tests.csv"))
count_labels <- table(model$model_group)
cplot <- ggplot(order_data, aes(Order, percent, fill = model_group)) +
  geom_boxplot(color = "black", position = position_dodge(.8), width = .65) +
  geom_segment(data = tests, aes(x = x - .2, xend = x + .2, y = 108, yend = 108), inherit.aes = FALSE) +
  geom_text(data = tests, aes(x = x, y = 113, label = label), inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = pal[c("DUSA-like", "Mesoplasma")],
                    labels = paste0(names(count_labels), " (n=", as.integer(count_labels), ")")) +
  coord_cartesian(ylim = c(0, 118)) + theme_classic(base_size = 12) +
  labs(x = NULL, y = "Relative read abundance (%)", fill = "Tested bacterial state",
       title = "Prey read abundance by order")
fig <- (a$plot + b$plot) / cplot + plot_layout(heights = c(1, 1.35)) + plot_annotation(
  tag_levels = "A")
ggsave(file.path(root, "figures/main/fig05_dominant_bacteria_nmds_boxplot.png"), fig,
       width = 16, height = 13, dpi = 300, bg = "white")
ggsave(file.path(root, "figures/main/fig05_dominant_bacteria_nmds_boxplot.pdf"), fig,
       width = 16, height = 13, device = "pdf", useDingbats = FALSE, bg = "white")
# S6 uses exactly the same 26 clear paired specimens and distance matrix as Figure5 tests.
dispersion_values <- read_csv(file.path(out, "16S_dispersion_values.csv"), show_col_types = FALSE)
s6 <- ggplot(dispersion_values, aes(group, distance_to_spatial_median, fill = group)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(shape = 21, color = "black", size = 2.5, stroke = .3,
              position = position_jitter(width = .1, seed = 109)) +
  scale_fill_manual(values = pal[c("DUSA-like", "Mesoplasma")]) +
  theme_classic() + theme(legend.position = "none") + labs(x = NULL, y = "Distance to spatial median")
ggsave(file.path(root, "figures/supplement/S06_betadispar.png"), s6,
       width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path(root, "figures/supplement/S06_betadispar.pdf"), s6,
       width = 6, height = 5, device = "pdf", useDingbats = FALSE, bg = "white")
write_csv(a$statistics %>% transmute(analysis = "microbiome_betadisper_main", n = n_model,
  n_DUSA = n_group1, n_Mesoplasma = n_group2, F = dispersion_F, p_value = dispersion_p,
  center = dispersion_center), file.path(root, "output/statistics/04_microbiome_betadisper_main.csv"))
write_csv(read_csv(file.path(out, "16S_ordination_scores.csv"), show_col_types = FALSE) %>%
  rename(arg_label = figure4_label), file.path(root, "output/tables/fig05_microbiome_nmds_scores.csv"))
write_csv(read_csv(file.path(out, "COI_ordination_scores.csv"), show_col_types = FALSE) %>%
  rename(arg_label = figure4_label), file.path(root, "output/tables/fig05_prey_nmds_scores.csv"))
write_csv(a$statistics, file.path(root, "output/statistics/fig05_microbiome_adonis.csv"))
write_csv(b$statistics, file.path(root, "output/statistics/fig05_prey_adonis.csv"))
write_csv(tests, file.path(root, "output/statistics/fig05_prey_order_WMW_BH.csv"))
capture.output(sessionInfo(), file = file.path(out, "session_info.txt"))
saveRDS(list(options = opt, command = commandArgs(), seeds = seeds, rng_kind = RNGkind(), final_rng_state = .Random.seed,
             cohort = cohort, rule = "Focal-only: Mixed when each of DUSA and Mesoplasma supplies at least 20% of their combined reads; otherwise the predominant focal group. Other bacterial reads do not determine grouping. Paired bacterial eligibility >=1500",
             statistical_input_order = "Immutable specimen ID, ascending; independent of display labels",
             source_files = source_files, input_md5 = tools::md5sum(source_files),
             iterations = 100L, permutations = 999L, figure_written = "figures/main/fig05_dominant_bacteria_nmds_boxplot.png"), file.path(out, "run_provenance.rds"))
cat("Figure 5 and S6: 28 paired profiles displayed; 26 clear profiles tested.\n")
