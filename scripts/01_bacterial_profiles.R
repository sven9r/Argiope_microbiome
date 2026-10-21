source(here::here("scripts", "00_setup.R"))

raw <- load_argiope_raw()
bacteria <- raw$bacteria
bacteria_taxa <- raw$bacteria_taxa

set.seed(100)
state <- make_dominance_state(bacteria, bacteria_taxa)
state_tbl <- state$state_tbl
cont_ab <- state$cont_ab
write_state_tables(state)

state_display <- function(x) {
  case_when(
    x == "ambiguous" ~ "Ambiguous",
    x == "DUSA_dominant" ~ "DUSA",
    x == "Mesoplasma_dominant" ~ "*Mesoplasma sp.*",
    TRUE ~ x
  )
}

fig04_dat <- make_fig04_arg_map(state_tbl, cont_ab) %>%
  mutate(
    other_rel = pmax(1 - dusa_rel - mesoplasma_rel, 0),
    dominance_display = state_display(dominance_state),
    repel_label = if_else(ID %in% c("TH81", "TH110"), arg_label, NA_character_)
  )

arg_levels <- fig04_dat$arg_label

s04_logratio_dat <- state_tbl %>%
  inner_join(
    fig04_dat %>% select(ID, arg_label),
    by = "ID"
  ) %>%
  inner_join(raw$prey %>% select(ID, PreyReads), by = "ID") %>%
  filter(PreyReads >= 1693) %>%
  arrange(match(arg_label, arg_levels)) %>%
  mutate(
    arg_label = factor(arg_label, levels = as.character(arg_label)),
    dominance_display = factor(state_display(dominance_state), levels = names(dominance_palette))
  )

stopifnot(nrow(s04_logratio_dat) == 28L, !any(s04_logratio_dat$ID %in% c("TH84", "TH99")))

readr::write_csv(
  s04_logratio_dat,
  here::here("output", "tables", "S04_logratio_lollipop_sample_order.csv")
)

s03_logratio <- s04_logratio_dat %>%
  ggplot(aes(arg_label, log10_dusa_mesoplasma, color = dominance_display)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_segment(aes(xend = arg_label, y = 0, yend = log10_dusa_mesoplasma), linewidth = 0.4) +
  geom_point(size = 2) +
  scale_color_manual(values = dominance_palette, breaks = names(dominance_palette), labels = c("Mixed", "DUSA", "*Mesoplasma*")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.text = element_markdown()) +
  labs(x = NULL, y = "log10[(DUSA-like + 1) / (Mesoplasma + 1)]", color = "Profile group")
save_supplement_figure(s03_logratio, "S05_logratio_lollipop.png", width = 12, height = 5)

meta_state <- state_tbl %>% filter(dominance_state %in% c("DUSA_dominant", "Mesoplasma_dominant", "ambiguous"))
bac_mat <- cont_ab %>%
  select(ID, all_of(state$bacteria_otu)) %>%
  as.data.frame()
rownames(bac_mat) <- bac_mat$ID
bac_mat <- as.matrix(bac_mat[, -1, drop = FALSE])

main_ids <- meta_state %>%
  filter(!is_TH81_or_TH110, dominance_state != "ambiguous") %>%
  pull(ID)
main_m <- bac_mat[intersect(main_ids, rownames(bac_mat)), , drop = FALSE]
main_meta <- meta_state %>% filter(ID %in% rownames(main_m))
main_m <- main_m[main_meta$ID, , drop = FALSE]
main_bc <- vegan::avgdist(main_m, sample = 1500, dmethod = "bray", iterations = 100)
labels_main <- attr(main_bc, "Labels")
if (!is.null(labels_main)) {
  main_meta <- main_meta %>% filter(ID %in% labels_main)
  main_meta <- main_meta[match(labels_main, main_meta$ID), , drop = FALSE]
}

bd <- vegan::betadisper(main_bc, main_meta$dominance_state, type = "median")
bd_p <- vegan::permutest(bd, permutations = 999)
readr::write_csv(
  tibble(
    analysis = "microbiome_betadisper_main",
    p_value = bd_p$tab[1, "Pr(>F)"]
  ),
  here::here("output", "statistics", "04_microbiome_betadisper_main.csv")
)

s04_betadisper <- tibble(distance_to_spatial_median = bd$distances, dominance_state = bd$group) %>%
  mutate(dominance_display = factor(state_display(dominance_state), levels = names(dominance_palette))) %>%
  ggplot(aes(dominance_display, distance_to_spatial_median, fill = dominance_display)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(shape = 21, width = 0.1, color = "black", size = 2.5, stroke = 0.3) +
  scale_fill_manual(values = dominance_palette, breaks = names(dominance_palette)) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_markdown()) +
  labs(x = NULL, y = "Distance to spatial median")
save_supplement_figure(s04_betadisper, "S06_betadispar.png", width = 6, height = 5)

log_line("01_bacterial_profiles completed")
