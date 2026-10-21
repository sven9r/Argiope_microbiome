source(here::here("scripts", "00_setup.R"))
raw <- load_argiope_raw()
state <- make_dominance_state(raw$bacteria, raw$bacteria_taxa)
write_state_tables(state)
cohort <- contemporary_analysis_cohort(state, raw$prey)
logratio <- cohort %>% filter(paired) %>%
  mutate(arg_label = factor(arg_label, levels = arg_label),
         profile_group = factor(recode(dominance_state, DUSA_dominant = "DUSA-like",
           Mesoplasma_dominant = "Mesoplasma", ambiguous = "Mixed"),
           levels = c("Mixed", "DUSA-like", "Mesoplasma")))
readr::write_csv(logratio, here::here("output/tables/S04_logratio_lollipop_sample_order.csv"))
s5 <- ggplot(logratio, aes(arg_label, log10_dusa_mesoplasma, color = profile_group)) +
  geom_hline(yintercept = 0, linewidth = .3, color = "grey40") +
  geom_segment(aes(xend = arg_label, y = 0, yend = log10_dusa_mesoplasma), linewidth = .4) +
  geom_point(size = 2) + scale_color_manual(values = c(Mixed = "black", `DUSA-like` = "#ff6a0a", Mesoplasma = "#eeee00")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
                         plot.margin = margin(8, 12, 20, 8)) +
  labs(x = NULL, y = "log10[(DUSA-like + 1) / (Mesoplasma + 1)]", color = "Profile group")
save_supplement_figure(s5, "S05_logratio_lollipop.png", width = 12, height = 5.4)
ggsave(here::here("figures/supplement/S05_logratio_lollipop.pdf"), s5,
       width = 12, height = 5.4, device = "pdf", useDingbats = FALSE, bg = "white")
log_line("S05_bacterial_profiles completed")
