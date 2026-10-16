thres_phylo <- (sqrt(max(all_SES_pair_long_clim$phy.cor))/3)^2
thres_phylo2 <- (sqrt(max(all_SES_pair_long_clim$phy.cor))/3*2)^2

example_df %>%
  mutate(
    native_sim_bin = case_when(
      native_sim < 0.1  ~ "Low (<0.1)",
      native_sim >= 0.4 & native_sim <= 0.6 ~ "Mid (0.4–0.6)",
      native_sim >= 0.9 ~ "High (>0.9)"
    ),
    native_sim_bin = factor(native_sim_bin, levels = c("Low (<0.1)", "Mid (0.4–0.6)", "High (>0.9)")),
    phy_cor_bin = case_when(
      phy.cor < thres_phylo  ~ "Low (<30 Mya)",
      phy.cor <= thres_phylo2 & phy.cor >= thres_phylo  ~ "Mid (>30 Mya and <120 Mya)",
      TRUE ~ "High (>120 Mya)"
    ),
    phy_cor_bin = factor(phy_cor_bin, levels = c("Low (<30 Mya)", "Mid (>30 Mya and <120 Mya)", "High (>120 Mya)"))
  ) %>%
  na.omit() %>%
  group_by(phy_cor_bin, native_sim_bin) %>%
  summarise(
    disp_y = sd(async),
    mean_y = mean(async),  
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    disp_y_scaled = (disp_y - min(disp_y)) / (max(disp_y) - min(disp_y)),
    text_color = ifelse(disp_y_scaled > 0.5, "black", "white")
  ) %>%
  ggplot(aes(x = phy_cor_bin, y = native_sim_bin, fill = disp_y)) +
  geom_tile() +
  geom_text(aes(label = paste0("Mean=", round(mean_y, 3), "\nSD=", round(disp_y, 3), "\nN=", n),
                color = text_color),
            size = 3.5) +
  scale_fill_viridis_c(name = "Interspecific asynchrony", option = "C") +
  scale_color_identity() +
  guides(fill = guide_colorbar(barheight = 1, barwidth = 8)) +
  theme_classic() +
  theme(legend.position="bottom") +
  labs(y = "Native distribution overlap", x = "Phylogenetic distance")

ggsave("Figure/Result/async_SI.tiff",
       dpi = 600, width = 6, height = 6, compression = "lzw")

example_df %>%
  mutate(
    native_sim_bin = case_when(
      native_sim < 0.1  ~ "Low (<0.1)",
      native_sim >= 0.4 & native_sim <= 0.6 ~ "Mid (0.4–0.6)",
      native_sim >= 0.9 ~ "High (>0.9)"
    ),
    native_sim_bin = factor(native_sim_bin, levels = c("Low (<0.1)", "Mid (0.4–0.6)", "High (>0.9)")),
    phy_cor_bin = case_when(
      phy.cor < thres_phylo  ~ "Low (<30 Mya)",
      phy.cor <= thres_phylo2 & phy.cor >= thres_phylo  ~ "Mid (>30 Mya and <120 Mya)",
      TRUE ~ "High (>120 Mya)"
    ),
    phy_cor_bin = factor(phy_cor_bin, levels = c("Low (<30 Mya)", "Mid (>30 Mya and <120 Mya)", "High (>120 Mya)"))
  ) %>%
  na.omit() %>%
  group_by(phy_cor_bin, native_sim_bin) %>%
  summarise(
    disp_y = sd(pop_stab_diff),
    mean_y = mean(pop_stab_diff),  
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    disp_y_scaled = (disp_y - min(disp_y)) / (max(disp_y) - min(disp_y)),
    text_color = ifelse(disp_y_scaled > 0.5, "black", "white")
  ) %>%
  ggplot(aes(x = phy_cor_bin, y = native_sim_bin, fill = disp_y)) +
  geom_tile() +
  geom_text(aes(label = paste0("Mean=", round(mean_y, 3), "\nSD=", round(disp_y, 3), "\nN=", n),
                color = text_color),
            size = 3.5) +
  scale_fill_viridis_c(name = "Population stability difference", option = "C") +
  scale_color_identity() +
  guides(fill = guide_colorbar(barheight = 1, barwidth = 8)) +
  theme_classic() +
  theme(legend.position="bottom") +
  labs(y = "Native distribution overlap", x = "Phylogenetic distance")


ggsave("Figure/Result/pop_stab_SI.tiff",
       dpi = 600, width = 6, height = 6, compression = "lzw")
