# Marina Natividad Avila
# 2024-12-06
# Final bar plot for GALA paper regarding neptune results for 73 acmg actionable genes

rm(list=ls())
`%notin%` <- Negate(`%in%`)

# Load necessary libraries
library(ggplot2)
library(ggsignif)

color_scheme <- c('#F3E9DC', '#9D6381', '#612940', '#C08552')
data <- data.frame(
  Comparison = rep(c('Neptune Classified/Total', 'PLP/Total', 'PLP/Neptune Classified'), each = 4),
  Group = rep(c("nonEUR", "AMR", "nonAMR", "EUR"), 3),
  PointEstimate = c(0.724, 0.699, 0.734, 0.723,   # Classified/Total
                    0.010, 0.0112, 0.0122, 0.014, # PLP/Total
                    0.013, 0.016, 0.0166, 0.019), # PLP/Neptune
  CI_lower = c(0.718, 0.691, 0.729, 0.717,       # Classified/Total CI Lower
               0.008, 0.009, 0.011, 0.013,       # PLP/Total CI Lower
               0.012, 0.013, 0.015, 0.017),      # PLP/Neptune CI Lower
  CI_upper = c(0.730, 0.707, 0.739, 0.729,       # Classified/Total CI Upper
               0.011, 0.013, 0.014, 0.016,       # PLP/Total CI Upper
               0.015, 0.019, 0.018, 0.022),      # PLP/Neptune CI Upper
  p_value = c(0.8153, 4.165E-13, 4.165E-13, 0.8153,   # Classified/Total
              4.87E-05, 0.4282, 0.4282, 4.87E-05,     # PLP/Total
              4.35E-05, 0.7612, 0.7612, 4.35E-05)     # PLP/Neptune
)

data$Group <- factor(data$Group, levels = c('nonEUR', 'AMR', 'nonAMR', 'EUR'))
data$Comparison <- factor(data$Comparison, levels = c('Neptune Classified/Total', 'PLP/Total', 'PLP/Neptune Classified'))

################################################################################
################################################################################
# Neptune Classified/Total
p1a <- ggplot(subset(data, Comparison == 'Neptune Classified/Total'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.7), size = 0.2) +
  labs(y = '', title = 'Neptune Classified/Total') +
  scale_y_continuous(breaks = c(seq(0, 0.90, 0.30)), limits = c(0, 0.90), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'), 
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('AMR', 'EUR'), c('AMR', 'nonAMR')),
              annotations = '*', y_position = c(0.84, 0.78), tip_length = 0, textsize = 1.5, size = 0.2, vjust = 0.3)

p1a

# PLP/Total
p1b <- ggplot(subset(data, Comparison == 'PLP/Total'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.7), size = 0.2) +
  labs(y = 'Point Estimate', title = 'PLP/Total') +
  scale_y_continuous(breaks = seq(0, 0.018, 0.006), limits = c(0, 0.018), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'), 
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('nonEUR', 'EUR'), c('AMR', 'EUR')),
              annotations = '*', y_position = c(0.017, 0.016), tip_length = 0, textsize = 1.5, size = 0.2, vjust = 0.3)

p1b

 # PLP/Neptune Classified
p1c <- ggplot(subset(data, Comparison == 'PLP/Neptune Classified'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.7), size = 0.2) +
  labs(y = '', title = 'PLP/Neptune Classified') +
  scale_y_continuous(breaks = seq(0, 0.027, 0.009), limits = c(0, 0.027), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'), 
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'),
        axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('nonEUR', 'EUR'), c('AMR', 'EUR')),
              annotations = '*', y_position = c(0.025, 0.023), tip_length = 0, textsize = 1.5, size = 0.2, vjust = 0.3)

p1c

fig7a <- cowplot::plot_grid(p1a, p1b, p1c, nrow = 3)
fig7a

################################################################################
################################################################################
################################################################################
data <- data.frame(
  Comparison = rep(c('Variant/Proband', 'Neptune Classified/Proband', 'PLP/Proband'), each = 4),
  Group = rep(c("AMR", "nonAMR", "nonEUR", "EUR"), 3),
  Total_Variants = c(12162, 28262, 19852, 20572,     # Proband/Total
                     8501, 20750, 14376, 14875,      # Proband/Neptune
                     136, 344, 191, 289),            # PLP/Proband
  N_Probands = c(4450, 13030, 6365, 11115,              # Proband/Total
              4450, 13030, 6365, 11115,              # Proband/Neptune
              4450, 13030, 6365, 11115),             # PLP/Proband
  PointEstimate = c(2.733, 2.169, 3.119, 1.851,      # Proband/Total
                    1.910, 1.592, 2.259, 1.338,      # Proband/Neptune
                    0.031, 0.026, 0.030, 0.026),     # PLP/Proband
  CI_lower = c(2.684, 2.144, 3.076, 1.826,           # Proband/Total CI Lower
               1.870, 1.571, 2.222, 1.317,           # Proband/Neptune CI Lower
               0.026, 0.024, 0.026, 0.023),          # PLP/Proband CI Lower
  CI_upper = c(2.782, 2.194, 3.162, 1.876,           # Proband/Total CI Upper
               1.951, 1.614, 2.296, 1.360,           # Proband/Neptune CI Upper
               0.036, 0.029, 0.035, 0.029),          # PLP/Proband CI Upper
  p_value = c(1.2549e-90, 1.2549e-90, 0, 0,          # Proband/Total p-values
              9.72e-42, 9.72e-42, 0, 0,              # Proband/Neptune p-values
              0.1575, 0.1575, 0.1306, 0.1306)        # PLP/Proband p-values
)

data$Group <- factor(data$Group, levels = c('nonEUR', 'AMR', 'nonAMR', 'EUR'))
data$Comparison <- factor(data$Comparison, levels = c('Variant/Proband', 'Neptune Classified/Proband', 'PLP/Proband'))

# Variant/Proband
p2a <- ggplot(subset(data, Comparison == 'Variant/Proband'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.8), size = 0.2) +
  labs(y = '', title = 'Variant/Proband') +
  scale_y_continuous(breaks = seq(0, 4.0, 1.0), limits = c(0, 4.0), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'),
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('nonEUR', 'EUR'), c('AMR', 'EUR'), c('AMR', 'nonAMR')),
              annotations = c('*', '*', '*'), y_position = c(3.6, 3.4, 3.2), tip_length = 0, textsize = 1.5, size = 0.2, vjust = 0.3)

p2a

# Neptune Classified/Proband
p2b <- ggplot(subset(data, Comparison == 'Neptune Classified/Proband'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.8), size = 0.2) +
  labs(y = '', title = 'Neptune Classified/Proband') +
  scale_y_continuous(breaks = seq(0, 2.8, 0.7), limits = c(0, 2.8), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'),
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('nonEUR', 'EUR'), c('AMR', 'EUR'), c('AMR', 'nonAMR')),
              annotations = c('*', '*', '*'), y_position = c(2.65, 2.5, 2.35), tip_length = 0, textsize = 1.5, size = 0.2, vjust = 0.3)

p2b

# PLP/Proband
p2c <- ggplot(subset(data, Comparison == 'PLP/Proband'), aes(x = Group, y = PointEstimate, fill = Group)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1, position = position_dodge(0.8), size = 0.2) +
  labs(y = '', title = 'PLP/Proband') +
  scale_y_continuous(breaks = seq(0, 0.042, 0.014), limits = c(0, 0.042), expand = c(0, 0)) +
  theme_classic() +
  scale_fill_manual(values = color_scheme) +
  theme(axis.text = element_text(size = 5, family = 'sans', color = 'black'), axis.title.y = element_text(size = 5, family = 'sans'),
        axis.title.x = element_blank(), title = element_text(size = 5, family = 'sans'), axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  geom_signif(comparisons = list(c('nonEUR', 'EUR'), c('AMR', 'EUR'), c('AMR', 'nonAMR')),
              annotations = c('p = 0.1306', 'p = 0.1329', 'p = 0.1575'), 
              y_position = c(0.040, 0.038, 0.036), tip_length = 0, textsize = 0.8, size = 0.2, vjust = 0.40)

p2c


fig7b <- cowplot::plot_grid(p2a, p2b, p2c, nrow = 3)
fig7b

fig7 <- cowplot::plot_grid(p1a, p2a, p1b, p2b, p1c, p2c, nrow = 3)
fig7

# ggsave('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Figures/manuscript/NatureMedicine/December/Supplementary_Figure4_2024-12-06.pdf', fig7, height = 120, width = 80, units = 'mm', device = cairo_pdf)

