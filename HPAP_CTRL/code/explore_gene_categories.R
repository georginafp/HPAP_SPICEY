# Plot distributions of number of enhancers per gene
reg_counts <- bind_rows(
  SPICEY_ANNOTATED_COACC %>%
    filter(genes_coacc %in% ct_markers$Gene) %>%
    mutate(type = "Tissue-specific"),

  SPICEY_ANNOTATED_COACC %>%
    filter(genes_coacc %in% ubiquitous) %>%
    mutate(type = "Ubiquitous")
) %>%
  mutate(type = factor(type, levels = c("Ubiquitous", "Tissue-specific"))) %>%
  distinct(genes_coacc, seqnames, start, end, type) %>%
  group_by(genes_coacc, type) %>%
  summarise(num_reg_elements = n(), .groups = "drop")

ymax <- max(reg_counts$num_reg_elements, na.rm = TRUE) * 1.05

ggplot(reg_counts, aes(x = type, y = num_reg_elements, fill = type)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = c("grey60", "steelblue")) +
  labs(x = NULL, y = "# Regulatory Elements") +
  theme_gray(base_size = 12) +
  theme(legend.position = "none") +
  stat_compare_means(
    method = "wilcox.test",
    label.x = 1.5,
    label.y = ymax
  )

# Test correlation between number of enhancers and entropy per gene
# Consider looking at strength of enhancers (RETSI values) weighted by entropy
# Check if some enhancers show bimodal patterns that affect entropy calculations


library(dplyr)
library(tidyr)
library(entropy)  # for entropy calculation

# Prepare data: remove NAs and keep relevant columns
df_filtered <- combined_df %>%
  filter(!is.na(genes_coacc), !is.na(RETSI), !is.na(GETSI_coacc)) %>%
  select(genes_coacc, cell_type, RETSI, GETSI = GETSI_coacc, type)

# Calculate entropy per gene and type for RETSI and GETSI
entropy_df <- df_filtered %>%
  group_by(genes_coacc, type) %>%
  summarise(
    RETSI_entropy = entropy::entropy(prop.table(RETSI + 1e-10)),  # add small offset to avoid zero probs
    GETSI_entropy = entropy::entropy(prop.table(GETSI + 1e-10)),
    n_enhancers = n_distinct(paste0(cell_type)),  # count distinct cell types (proxy for enhancer diversity)
    .groups = "drop"
  )


num_re_df <- combined_df %>%
  distinct(genes_coacc, seqnames, start, end, type) %>%
  group_by(genes_coacc, type) %>%
  summarise(num_reg_elements = n(), .groups = "drop")


merged_df <- entropy_df %>%
  inner_join(num_re_df, by = c("genes_coacc", "type"))



cor_retsi <- cor.test(merged_df$num_reg_elements, merged_df$RETSI_entropy, method = "spearman")
cor_getsi <- cor.test(merged_df$num_reg_elements, merged_df$GETSI_entropy, method = "spearman")


p1 <- ggplot(merged_df, aes(x = num_reg_elements, y = RETSI_entropy, color = type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Correlation between number of enhancers and RETSI entropy",
       x = "# Regulatory Elements per Gene",
       y = "RETSI Entropy") +
  theme_minimal()

p2 <- ggplot(merged_df, aes(x = num_reg_elements, y = GETSI_entropy, color = type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Correlation between number of enhancers and GETSI entropy",
       x = "# Regulatory Elements per Gene",
       y = "GETSI Entropy") +
  theme_minimal()

library(gridExtra)
gridExtra::grid.arrange(p1, p2, ncol = 2)



ggplot(merged_df, aes(x = type, y = RETSI_entropy, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "RETSI Entropy by Gene Category", y = "RETSI Entropy") +
  scale_fill_manual(values = c("grey60", "steelblue")) +
  theme_minimal()

ggplot(merged_df, aes(x = type, y = GETSI_entropy, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "GETSI Entropy by Gene Category", y = "GETSI Entropy") +
  scale_fill_manual(values = c("grey60", "steelblue")) +
  theme_minimal()

