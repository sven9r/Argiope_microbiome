source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
prey <- raw$prey
prey_taxa <- raw$prey_taxa
bacteria <- raw$bacteria
bacteria_taxa <- raw$bacteria_taxa

set.seed(100)
state <- make_dominance_state(bacteria, bacteria_taxa)
state_tbl <- state$state_tbl
fig04_arg_map <- make_fig04_arg_map(state_tbl, state$cont_ab)

coi_otu <- grep("^OTU", names(prey), value = TRUE)

fig05_prey_samples <- prey %>%
  filter(
    Status == "Contemporary",
    str_detect(str_to_lower(Taxon), "bruennichi"),
    PreyReads >= 1693
  ) %>%
  left_join(state_tbl %>% select(ID, dominance_state, is_TH81_or_TH110), by = "ID") %>%
  arrange(DominantSymbiont, Number) %>%
  apply_fig04_arg_labels(fig04_arg_map, fallback = "next_arg")

s06_sample_order <- fig05_prey_samples %>%
  mutate(fig04_order = match(arg_label, fig04_arg_map$arg_label)) %>%
  arrange(is.na(fig04_order), fig04_order, DominantSymbiont, Number) %>%
  transmute(ID, spider_label = arg_label)

readr::write_csv(
  fig05_prey_samples %>% select(ID, Number, Taxon, PreyReads, DominantSymbiont, dominance_state),
  here::here("output", "tables", "S07_detailed_diet_sample_set.csv")
)

diet_orders <- fig05_prey_samples %>%
  select(ID, Number, DominantSymbiont, all_of(coi_otu)) %>%
  pivot_longer(cols = all_of(coi_otu), names_to = "OTU", values_to = "reads") %>%
  left_join(prey_taxa %>% select(OTU, Order), by = "OTU") %>%
  mutate(Order = if_else(is.na(Order) | Order == "", "Unknown", Order)) %>%
  group_by(ID, Number, DominantSymbiont, Order) %>%
  summarize(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
  group_by(ID) %>%
  mutate(rel_abund = 100 * reads / pmax(sum(reads, na.rm = TRUE), 1)) %>%
  ungroup() %>%
  mutate(
    Order = if_else(
      Order %in% c("Blattodea", "Coleoptera", "Diptera", "Hemiptera", "Hymenoptera", "Lepidoptera", "Orthoptera"),
      Order,
      "Other (<3%)"
    )
  ) %>%
  group_by(ID, Number, DominantSymbiont, Order) %>%
  summarize(rel_abund = sum(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  left_join(s06_sample_order, by = "ID") %>%
  arrange(DominantSymbiont, Number) %>%
  mutate(
    spider_label = factor(spider_label, levels = unique(spider_label)),
    Order = factor(Order, levels = names(prey_order_palette)[names(prey_order_palette) %in% unique(Order)])
  )

readr::write_csv(
  diet_orders,
  here::here("output", "tables", "S07_detailed_diet_stack_data.csv")
)

matched_ids <- intersect(fig05_prey_samples$ID, state$cont_ab$ID)
bac_m <- state$cont_ab %>%
  filter(ID %in% matched_ids) %>%
  select(ID, all_of(state$bacteria_otu)) %>%
  as.data.frame()
rownames(bac_m) <- bac_m$ID
bac_m <- as.matrix(bac_m[, -1, drop = FALSE])

prey_m <- fig05_prey_samples %>%
  filter(ID %in% matched_ids) %>%
  select(ID, all_of(coi_otu)) %>%
  as.data.frame()
rownames(prey_m) <- prey_m$ID
prey_m <- as.matrix(prey_m[, -1, drop = FALSE])

common_ids <- sort(intersect(rownames(bac_m), rownames(prey_m)))
bac_m <- bac_m[common_ids, , drop = FALSE]
prey_m <- prey_m[common_ids, , drop = FALSE]

bac_dist <- vegan::avgdist(bac_m, sample = 1500, dmethod = "bray", iterations = 100)
prey_dist <- vegan::avgdist(prey_m, sample = 1693, dmethod = "bray", iterations = 100)

bac_labels <- attr(bac_dist, "Labels")
prey_labels <- attr(prey_dist, "Labels")
if (!is.null(bac_labels) && !is.null(prey_labels)) {
  keep <- intersect(bac_labels, prey_labels)
  bac_dist <- as.dist(as.matrix(bac_dist)[keep, keep, drop = FALSE])
  prey_dist <- as.dist(as.matrix(prey_dist)[keep, keep, drop = FALSE])
}

mantel_res <- vegan::mantel(bac_dist, prey_dist, method = "spearman", permutations = 999)
readr::write_csv(
  tibble(
    version = "fig05_sample_set_abundance",
    n = attr(bac_dist, "Size"),
    mantel_r = unname(mantel_res$statistic),
    p_value = mantel_res$signif
  ),
  here::here("output", "statistics", "S09_mantel_summary.csv")
)

bac_mat <- as.matrix(bac_dist)
prey_mat <- as.matrix(prey_dist)
keep <- intersect(rownames(bac_mat), rownames(prey_mat))
bac_mat <- bac_mat[keep, keep, drop = FALSE]
prey_mat <- prey_mat[keep, keep, drop = FALSE]
idx <- upper.tri(bac_mat)
mantel_plot_dat <- tibble(
  prey_dissimilarity = prey_mat[idx],
  microbiome_dissimilarity = bac_mat[idx]
)

s05_mantel <- ggplot(mantel_plot_dat, aes(prey_dissimilarity, microbiome_dissimilarity)) +
  geom_point(alpha = 0.55, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black") +
  theme_classic() +
  labs(
    x = "Prey Bray-Curtis dissimilarity",
    y = "Microbiome Bray-Curtis dissimilarity",
    title = sprintf(
      "Mantel: rho = %.3f, p = %.3f",
      unname(mantel_res$statistic),
      mantel_res$signif
    )
  )

save_supplement_figure(s05_mantel, "S09_mantel.png", width = 7, height = 6)

log_line("02_prey_profiles completed")
