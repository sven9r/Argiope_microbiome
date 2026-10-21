#!/usr/bin/env Rscript
# Supplementary Figure 4: species-level bacterial and prey read composition.
source(here::here("scripts", "00_setup.R"))
raw <- load_argiope_raw()
bacteria <- filter(raw$bacteria, !is.na(ID))
prey <- filter(raw$prey, !is.na(ID))

count_matrix <- function(data, read_column) {
  matrix <- as.matrix(select(data, starts_with("OTU")))
  storage.mode(matrix) <- "integer"
  rownames(matrix) <- data$ID
  stopifnot(!anyDuplicated(data$ID), !anyNA(matrix), all(matrix >= 0),
            all(rowSums(matrix) == data[[read_column]]))
  matrix
}
bacterial_counts <- count_matrix(bacteria, "BacteriaReads")
prey_counts <- count_matrix(prey, "PreyReads")

# The original analysis rarefies each marker once, resetting seed 555.
# Museum and contemporary specimens are pooled within host species.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(555)
bacterial_rarefied <- GUniFrac::Rarefy(bacterial_counts, depth = 507)$otu.tab.rff
set.seed(555)
prey_rarefied <- GUniFrac::Rarefy(prey_counts, depth = 1133)$otu.tab.rff
stopifnot(all(rowSums(bacterial_rarefied) == 507),
          all(rowSums(prey_rarefied) == 1133),
          setequal(rownames(bacterial_rarefied), bacteria$ID[bacteria$BacteriaReads >= 507]),
          setequal(rownames(prey_rarefied), prey$ID[prey$PreyReads >= 1133]))

# S4B additionally requires the workbook's original valid-prey classification.
bacterial_ids <- rownames(bacterial_rarefied)
prey_ids <- prey$ID[prey$PreyReads >= 1133 & prey$Category == "valid"]
prey_rarefied <- prey_rarefied[prey_ids, , drop = FALSE]

long_counts <- function(matrix, metadata) {
  as.data.frame(matrix) %>%
    rownames_to_column("ID") %>%
    pivot_longer(-ID, names_to = "OTU", values_to = "rarefied_reads") %>%
    left_join(select(metadata, ID, Taxon), by = "ID")
}
bacterial_long <- long_counts(bacterial_rarefied, bacteria) %>%
  left_join(select(raw$bacteria_taxa, OTU, Phylum), by = "OTU") %>%
  mutate(group = case_when(
    is.na(Phylum) | Phylum == "" ~ "Uncultured bacterium",
    Phylum == "Unknown symbiont" ~ "DUSA",
    TRUE ~ Phylum
  ))
prey_orders <- c("Coleoptera", "Diptera", "Hemiptera", "Hymenoptera",
                 "Lepidoptera", "Orthoptera")
prey_long <- long_counts(prey_rarefied, prey) %>%
  left_join(select(raw$prey_taxa, OTU, Order), by = "OTU") %>%
  mutate(group = if_else(!is.na(Order) & Order %in% prey_orders, Order, "Other"))

# Pool rarefied counts by species, then divide by that species' total reads.
species_composition <- function(data, marker, depth) {
  data %>%
    group_by(Taxon, group) %>%
    summarize(rarefied_reads = sum(rarefied_reads), .groups = "drop") %>%
    group_by(Taxon) %>%
    mutate(pooled_rarefied_reads = sum(rarefied_reads),
           relative_read_percent = 100 * rarefied_reads / pooled_rarefied_reads) %>%
    ungroup() %>%
    filter(rarefied_reads > 0) %>%
    mutate(marker = marker, rarefaction_depth = depth)
}
bacterial_stack <- species_composition(bacterial_long, "16S", 507)
prey_stack <- species_composition(prey_long, "COI", 1133)
bacterial_n <- bacteria %>% filter(ID %in% bacterial_ids) %>% count(Taxon, name = "n")
prey_n <- prey %>% filter(ID %in% prey_ids) %>% count(Taxon, name = "n")

stopifnot(sum(bacterial_n$n) == 56L, nrow(bacterial_n) == 19L,
          sum(prey_n$n) == 48L, nrow(prey_n) == 12L,
          bacterial_n$n[bacterial_n$Taxon == "A. bruennichi"] == 34L,
          prey_n$n[prey_n$Taxon == "A. bruennichi"] == 35L)
plot_data <- bind_rows(bacterial_stack, prey_stack)
totals <- plot_data %>% group_by(marker, Taxon) %>%
  summarize(total = sum(relative_read_percent), .groups = "drop")
stopifnot(all(abs(totals$total - 100) < 1e-10))

# Named colors retain the figure's category order and color assignments.
bacterial_colors <- c(
  Actinobacteria = "#99000D", Armatimonadetes = "#CB181D",
  Bacillariophyta = "#737373", Bacteroidota = "#FB6A4A",
  Candidatus = "#FF7F00", Chlamydiae = "#FFA500", Fibrobacteres = "#EEEE00",
  Firmicutes = "#74C476", Fusobacteria = "#238B45", Planctomycetota = "#005A32",
  Proteobacteria = "#6BAED6", Synergistetes = "#4292C6", Tenericutes = "#2171B5",
  DUSA = "#FF34B3", `Uncultured bacterium` = "#999999"
)
prey_colors <- c(
  Coleoptera = "#99000D", Diptera = "#FB6A4A", Hemiptera = "#EEEE00",
  Hymenoptera = "#74C476", Lepidoptera = "#6BAED6", Orthoptera = "#2171B5",
  Other = "#999999"
)
stopifnot(all(bacterial_stack$group %in% names(bacterial_colors)),
          all(prey_stack$group %in% names(prey_colors)))

composition_panel <- function(data, counts, colors, title, legend_title) {
  labels <- setNames(paste0(counts$Taxon, "\nn=", counts$n), counts$Taxon)
  data <- data %>%
    mutate(Taxon = factor(Taxon, levels = sort(unique(Taxon))),
           group = factor(group, levels = names(colors)))
  ggplot(data, aes(Taxon, relative_read_percent, fill = group)) +
    geom_col(width = .7, colour = "black") +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(title = title, x = NULL, y = "Relative read abundance (%)", fill = legend_title) +
    theme_classic() +
    theme(legend.text = element_text(size = 13),
          legend.key.size = grid::unit(12, "pt"),
          legend.title = element_text(size = 15, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 20),
          legend.spacing.y = grid::unit(.2, "cm")) +
    guides(fill = guide_legend(byrow = TRUE))
}
panel_a <- composition_panel(bacterial_stack, bacterial_n, bacterial_colors,
  "Relative abundance of the microbial gut communities\n", "Phylum")
panel_b <- composition_panel(prey_stack, prey_n, prey_colors,
  "Relative abundance of the main prey insect orders\n", "Order")
figure <- (panel_a / panel_b) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))
save_supplement_figure(figure, "S04_taxon_and_prey_stack.png", width = 14, height = 12)
write_csv(plot_data, here::here("output", "tables", "S04_species_composition.csv"))
write_csv(bind_rows(mutate(bacterial_n, marker = "16S"), mutate(prey_n, marker = "COI")),
          here::here("output", "tables", "S04_species_sample_counts.csv"))
cat(here::here("figures", "supplement", "S04_taxon_and_prey_stack.png"), "\n")
