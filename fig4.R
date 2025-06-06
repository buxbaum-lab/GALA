# Marina Natividad Avila
# 2024-10-03
# Figure 4 for Nature Genetics

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

rm(list=ls())
`%notin%` <- Negate(`%in%`)

metrics <- fread("/sc/arion/projects/buxbaj01a/GALA/DATA/final_data/gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt")

# TADA results
df <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/manuscript/2024/TADA_results/amr_bf_table_2024-04-17.txt')
# df <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/TADA-results/amr_bf_table_dn_inh_cc_2024-03-28.txt')
#df <- df[df$qval <= 0.01] # 25 on Oct 16

df$chr <- NA
df$chr <- metrics$chromosome[match(df$gene_id, metrics$gene_id)]
  df$chr <- as.numeric(df$chr)
df$gene_length <- NA
df$gene_length <- metrics$gene_length[match(df$gene_id, metrics$gene_id)]
df$start_position <- NA
df$start_position <- metrics$start_position[match(df$gene_id, metrics$gene_id)]
  df$start_position <- as.numeric(df$start_position)

df <- df[with(df, order(chr, start_position))]
df$idx <- 1:nrow(df)

don <- df %>% 
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(start_position)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, start_position) %>%
  mutate( BPcum=start_position+tot)


axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

manhattan.plt <- ggplot(don, aes(x = BPcum, y = -log10(qval_cc))) +
  # Show all points, scale down size
  geom_jitter(aes(color = as.factor(chr)), alpha = 0.8, size = 2 * 0.4375) +  # Scaled down points
  scale_x_continuous(label = axisdf$chr, breaks = axisdf$center, expand = c(0.005, 0.005)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 11, 2), limits = c(0, 11)) +
  theme_minimal() +
  #geom_hline(yintercept = -log10(0.1), linetype = 2, colour = '#CECCCC') +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = '#837E80') +
  geom_point(data = don[don$qval_cc < 0.01], color = "#612940", size = 4.5 * 0.4375) +  # Scaled down points
  geom_point(data = don[don$qval_cc < 0.05 & don$qval_cc >= 0.01], color = "#9D6381", size = 3 * 0.4375) +  # Scaled down points
  ggrepel::geom_text_repel(data = don[don$qval_cc < 0.01], aes(label = gene), size = 7 * 0.4, 
                           max.overlaps = Inf, fontface = 'italic', family = 'sans', box.padding = unit(0.4, 'mm')) +
  ggrepel::geom_text_repel(data = don[don$qval_cc < 0.05 & don$qval_cc >= 0.01], aes(label = gene), size = 4 * 0.4,
                           max.overlaps = Inf, fontface = 'italic', family = 'sans', box.padding = unit(0.4, 'mm')) +
  # ggrepel::geom_text_repel(data = don[don$qval_cc >= 0.01 & don$qval_cc < 0.1], aes(label = gene), size = 2.5 * 0.4, 
  #                          max.overlaps = Inf, fontface = 'italic', family = 'sans', segment.colour = NA, box.padding = unit(0.40, 'mm')) +
  labs(x = 'Chromosome', y = expression(paste(-log[10], '(Q)', sep = ''))) +
  scale_color_manual(values = rep(c('#CECCCC', '#D8A5AF'), 22)) +
  theme(
    axis.text = element_text(size = 7, family = 'sans', color = 'black'),
    axis.title = element_text(size = 7, family = 'sans'),
    plot.caption.position = 'plot', 
    plot.caption = element_text(hjust = 0, family = 'sans', size = 7),
    legend.position = 'none', 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    axis.ticks.length = unit(.15, 'cm'),
    axis.line = element_line(), 
    plot.margin = ggplot2::margin(0.25, 0.25, 0.25, 0.25, 'cm'),
    axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), family = 'sans', size = 7),
    axis.title.x = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), family = 'sans', size = 7),
    plot.title = element_text(size = 7, family = 'sans', hjust = 0.95),
    axis.ticks = element_line(lineend = 'square'))

manhattan.plt

# ggsave('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Figures/manuscript/NatureMedicine/December/Figure4_2024-12-06.pdf', manhattan.plt, height = 100, width = 180, units = 'mm')

################################################################################
################################################################################