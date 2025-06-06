# Marina Natividad Avila
# 2025-05-08
# redoing Figure 2 and Supplementary Figure 2

rm(list=ls())
`%notin%` <- Negate(`%in%`)

library(grid)
library(ggplot2)
library(data.table)

## ALL, AMR and notAMR
denovos <- openxlsx::read.xlsx('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/submissions/Supplementary_Tables.xlsx', sheet = 'Supplementary Table 2')
denovos <- setDT(denovos)

#################################################
#################################################

denovos[isPTV == TRUE & LOEUF_bin <= 2, ':=' (Tag2 = 'PTV.0')]
denovos[isPTV == TRUE & LOEUF_bin >= 3, ':=' (Tag2 = 'PTV.1')]
denovos[(isMis == TRUE & MPC >= 2), ':=' (Tag2 = 'MisB')]
denovos[(isMis == TRUE & MPC >= 1 & MPC < 2), ':=' (Tag2 = 'MisA')]
denovos[(isMis == TRUE & (MPC < 1 | is.na(MPC))), ':=' (Tag2 = 'MisX')]
denovos[isSyn == TRUE, ':=' (Tag2 = 'SYN')]

table(denovos$Tag, denovos$Tag2)

denovos$Pop2 <- denovos$Pop
denovos$Pop2[denovos$Pop2 != 'AMR'] <- 'OTH'

denovos.table <- as.data.frame(table(denovos$Role, denovos$Pop2, denovos$Tag2))
names(denovos.table) <- c('Affected_Status', 'Pop', 'Variant_Type', 'denovos')

all.table <- as.data.frame(table(denovos$Role, denovos$Tag2))
names(all.table) <- c('Affected_Status', 'Variant_Type', 'denovos')
all.table$Pop <- 'ALL'

all.table <- all.table[, c('Affected_Status', 'Pop', 'Variant_Type', 'denovos')]

fig.df <- rbind(all.table, denovos.table)

###############
### sample size
pedigree <- openxlsx::read.xlsx('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/submissions/Supplementary_Tables.xlsx', sheet = 'Supplementary Table 1')
with(pedigree, table(Affected_Status, Pop))
#                 Pop
# Affected_Status   AFR   AMR   EAS   FIN   MID   NFE   SAS
#               1   184  1459   112     7   121  4079   246
#               2   477  4450   419   117   446 10998   573

pedigree$Pop2 <- pedigree$Pop
pedigree$Pop2[pedigree$Pop2 != 'AMR'] <- 'OTH'
with(pedigree, table(Affected_Status, Pop2))
#               Pop2
# Affected_Status   AMR   OTH
#               1  1459  4749
#               2  4450 13030

sample.size <- as.data.frame(table(pedigree$Affected_Status, pedigree$Pop2))
names(sample.size) <- c('Affected_Status', 'Pop', 'N')

sample.size.all <- as.data.frame(table(pedigree$Affected_Status))
names(sample.size.all) <- c('Affected_Status', 'N')
sample.size.all$Pop <- 'ALL'
sample.size.all <- sample.size.all[, c('Affected_Status', 'Pop', 'N')]

sample.size.final <- rbind(sample.size, sample.size.all)

sample.size.final$Role <- NA
sample.size.final$Role[sample.size.final$Affected_Status == 1] <- 'Sibling'
sample.size.final$Role[sample.size.final$Affected_Status == 2] <- 'Proband'

table(sample.size.final$Role, sample.size.final$Affected_Status)

sample.size.final$key <- paste(sample.size.final$Role, sample.size.final$Pop, sep = '::')
###############
###############
fig.df$key <- paste(fig.df$Affected_Status, fig.df$Pop, sep = '::')

fig.df$N <- NA
fig.df$N <- sample.size.final$N[match(fig.df$key, sample.size.final$key)]
dim(fig.df[is.na(fig.df$N),])

fig.df$Rate <- fig.df$denovos/fig.df$N  

fig.df$CI.low <- qchisq(0.025, 2*fig.df$denovos)/2
fig.df$CI.high <- qchisq(1 - 0.025, 2*(fig.df$denovos + 1))/2

fig.df$CI.low.norm <- fig.df$CI.low/fig.df$N
fig.df$CI.high.norm <- fig.df$CI.high/fig.df$N

################################################################################
################################################################################
fig.df <- fig.df %>% group_by(Pop, Affected_Status) %>%
  mutate(SYN_rate = Rate[Variant_Type == 'SYN'],
         Rate_normalized = Rate / SYN_rate) %>% ungroup()

fig.df <- fig.df %>%
  group_by(Pop, Affected_Status) %>%
  mutate(
    SE = (CI.high.norm - Rate) / 1.96,
    SE_SYN = SE[Variant_Type == 'SYN'],
    Rate_SYN = Rate[Variant_Type == 'SYN'],
    Rate_normalized = Rate / Rate_SYN,
    SE_norm = sqrt((SE / Rate)^2 + (SE_SYN / Rate_SYN)^2) * Rate_normalized,
    CI.low.norm.adj = Rate_normalized - 1.96 * SE_norm,
    CI.high.norm.adj = Rate_normalized + 1.96 * SE_norm
  ) %>% ungroup()

################################################################################
################################################################################
# fig.df$Pop <- factor(fig.df$Pop, levels = c('ALL', 'AMR', 'NFE', 'AFR', 'SAS', 'EAS', 'MID', 'FIN'))
# fig.df$Variant_Type <- factor(fig.df$Variant_Type, levels = c('PTV.0', 'PTV.1', 'MisB', 'MisA', 'MisX', 'SYN'))

#fig.df$Rate.Log <- log(fig.df$Rate)

fig.df$variant_type <- as.character(fig.df$Variant_Type)

################################################################################
################################################################################
normalized_stats <- fig.df %>%
  filter(Affected_Status %in% c('Proband', 'Sibling')) %>%
  select(Pop, variant_type, Affected_Status, Rate_normalized, SE_norm) %>%
  pivot_wider(
    names_from = Affected_Status,
    values_from = c(Rate_normalized, SE_norm),
    names_sep = '.'
  ) %>%
  mutate(
    z = (Rate_normalized.Proband - Rate_normalized.Sibling) /
      sqrt(SE_norm.Proband^2 + SE_norm.Sibling^2),
    p_value = 2 * pnorm(-abs(z)),
    p_adjusted = p.adjust(p_value, method = 'fdr'),
    signif = ifelse(p_adjusted < 0.05, '*', 'ns')
  )

# View AMR results
normalized_stats
#   Pop   variant_type Rate_normalized.Proband Rate_normalized.Sibling SE_norm.Proband SE_norm.Sibling       z  p_value p_adjusted signif
# <chr> <chr>                          <dbl>                   <dbl>           <dbl>           <dbl>   <dbl>    <dbl>      <dbl> <chr> 
# 1 ALL   MisA                           0.488                   0.485         0.0123          0.0216   0.122  9.03e- 1   1   e+ 0 ns    
# 2 ALL   MisB                           0.200                   0.121         0.00718         0.00970  6.52   6.88e-11   4.13e-10 *     
# 3 ALL   MisX                           1.90                    1.93          0.0338          0.0599  -0.415  6.78e- 1   1   e+ 0 ns    
# 4 ALL   PTV.0                          0.295                   0.139         0.00901         0.0104  11.4    6.65e-30   1.20e-28 *     
# 5 ALL   PTV.1                          0.195                   0.196         0.00707         0.0125  -0.0281 9.78e- 1   1   e+ 0 ns    
# 6 ALL   SYN                            1                       1             0.0204          0.0357   0      1   e+ 0   1   e+ 0 ns    
# 7 AMR   MisA                           0.462                   0.527         0.0239          0.0496  -1.18   2.37e- 1   6.10e- 1 ns    
# 8 OTH   MisA                           0.497                   0.473         0.0145          0.0242   0.861  3.89e- 1   7.53e- 1 ns    
# 9 AMR   MisB                           0.198                   0.137         0.0145          0.0233   2.23   2.55e- 2   7.66e- 2 ns    
# 10 OTH   MisB                           0.201                   0.117         0.00838         0.0109   6.10   1.04e- 9   4.70e- 9 *     
# 11 AMR   MisX                           1.91                    2.03          0.0679          0.135   -0.809  4.18e- 1   7.53e- 1 ns    
# 12 OTH   MisX                           1.90                    1.90          0.0393          0.0674  -0.0238 9.81e- 1   1   e+ 0 ns    
# 13 AMR   PTV.0                          0.253                   0.142         0.0166          0.0237   3.81   1.41e- 4   5.09e- 4 *     
# 14 OTH   PTV.0                          0.310                   0.138         0.0108          0.0119  10.7    7.91e-27   7.12e-26 *     
# 15 AMR   PTV.1                          0.197                   0.228         0.0144          0.0303  -0.925  3.55e- 1   7.53e- 1 ns    
# 16 OTH   PTV.1                          0.194                   0.186         0.00822         0.0139   0.519  6.04e- 1   9.88e- 1 ns    
# 17 AMR   SYN                            1                       1             0.0408          0.0772   0      1   e+ 0   1   e+ 0 ns    
# 18 OTH   SYN                            1                       1             0.0237          0.0407   0      1   e+ 0   1   e+ 0 ns 
################################################################################
################################################################################

### previous results
# Pop   variant_type  p_value p_adjusted signif
# <chr> <chr>           <dbl>      <dbl> <chr> 
# 1 ALL   MisA         5.57e- 2   7.71e- 2 ns    
# 2 ALL   MisB         3.71e-14   2.23e-13 *     
# 3 ALL   MisX         2.09e- 4   6.26e- 4 *     
# 4 ALL   PTV.0        9.17e-35   1.65e-33 *     
# 5 ALL   PTV.1        3.14e- 1   3.53e- 1 ns    
# 6 ALL   SYN          5.70e- 3   1.14e- 2 *     
# 7 AMR   MisA         1.00e+ 0   1.00e+ 0 ns    
# 8 AMR   MisB         9.91e- 4   2.55e- 3 *     
# 9 AMR   MisX         1.93e- 2   3.16e- 2 *     
# 10 AMR   PTV.0        8.42e- 7   3.03e- 6 *     
# 11 AMR   PTV.1        9.38e- 1   9.93e- 1 ns    
# 12 AMR   SYN          1.11e- 2   2.00e- 2 *     
# 13 OTH   MisA         2.62e- 2   3.92e- 2 *     
# 14 OTH   MisB         1.12e-11   5.06e-11 *     
# 15 OTH   MisX         4.01e- 3   9.02e- 3 *     
# 16 OTH   PTV.0        1.06e-29   9.58e-29 *     
# 17 OTH   PTV.1        2.20e- 1   2.65e- 1 ns    
# 18 OTH   SYN          8.61e- 2   1.11e- 1 ns  
################################################################################
################################################################################
fig.df$Variant_Type <- factor(fig.df$Variant_Type, levels = c('PTV.0', 'PTV.1', 'MisB', 'MisA', 'MisX', 'SYN'))

################################################################################
################################################################################

## LOEUF
fig1a <- ggplot(fig.df[fig.df$Variant_Type %in% c('PTV.0', 'PTV.1'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('PTV.0' = 'PTV 1\u02e2\u1d57-3\u02b3\u1d48\n LOEUF deciles',
                                                                 'PTV.1' = 'PTV 4\u1d57\u02b0-10\u1d57\u02b0\n LOEUF deciles')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = '', labels = c('ALL' = 'ALL', 'AMR' = 'AMR', 'OTH' = expression(Fu[COMP]))) + ylab('Rate') +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 6, family = 'sans'), strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0.2, 0.2, -1, 0.2, 'cm'), panel.spacing = unit(0.1, 'lines'), legend.position = c(0.92, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 2, l = 0), size = 5, family = 'sans'), axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 20, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'PTV.0' ~ scale_y_continuous(breaks = c(seq(0, 0.36, 0.12)), limits = c(0, 0.36), expand = c(0, 0)),
    Variant_Type == 'PTV.1' ~ scale_y_continuous(breaks = c(seq(0, 0.36, 0.12)), limits = c(0, 0.36), expand = c(0, 0)))) +
  geom_signif(y_position = c(0.32, 0.23), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p < 1e-28', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(0.295, 0.32), xmin = c(1.8, 1.8), xmax = c(2.2, 2.2), annotation = c('p < 5e-4', 'p = 0.753'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(0.34, 0.22), xmin = c(2.8, 2.8), xmax = c(3.2, 3.2), annotation = c('p < 1e-26', 'p = 0.988'), tip_length = 0, textsize = 0.9, size = 0.2)

fig1a

 
## missense
fig1b <- ggplot(fig.df[fig.df$Variant_Type %in% c('MisB', 'MisA'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('MisB' = 'missense\n MPC \u2265 2',
                                                                 'MisA' = 'missense\n 1 \u2264 MPC < 2')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = 'Population', labels = c('ALL' = 'ALL', 'AMR' = 'AMR', 'OTH' = expression(Fu[COMP]))) +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 6, family = 'sans'), strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0.2, 0.2, -1, 0.2, 'cm'), panel.spacing = unit(0.1, 'lines'), legend.position = c(0.92, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 2, l = 0), size = 5, family = 'sans'), axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 20, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'MisB' ~ scale_y_continuous(breaks = c(seq(0, 0.70, 0.14)), limits = c(0, 0.70), name = '', expand = c(0, 0)),
    Variant_Type == 'MisA' ~ scale_y_continuous(breaks = c(seq(0, 0.70, 0.14)), limits = c(0, 0.70), name = '', expand = c(0, 0)))) +
  geom_signif(y_position = c(0.23, 0.55), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p < 1e-10', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(0.25, 0.64), xmin = c(1.8, 1.8), xmax = c(2.2, 2.2), annotation = c('p = 0.077', 'p = 0.61'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(0.24, 0.55), xmin = c(2.8, 2.8), xmax = c(3.2, 3.2), annotation = c('p < 1e-9', 'p = 0.753'), tip_length = 0, textsize = 0.9, size = 0.2)


fig1b

## synonymous
fig1c <- ggplot(fig.df[fig.df$Variant_Type %in% c('MisX', 'SYN'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = TRUE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('MisX' = 'missense\n MPC < 1',
                                                                 'SYN' = 'synonymous\n')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = '', labels = c('ALL' = 'ALL', 'AMR' = 'AMR', 'OTH' = expression(Fu[COMP]))) +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        legend.key.size = unit(0.2, 'cm'), legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0.1, 'cm'), 
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 6, family = 'sans'), strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0.2, 0.2, -1, 0.2, 'cm'), legend.position = c(0.85, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 2, l = 0), size = 5, family = 'sans'), axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 20, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'MisX' ~ scale_y_continuous(breaks = c(seq(0, 2.4, 0.60)), limits = c(0, 2.4), name = '', expand = c(0, 0)),
    Variant_Type == 'SYN' ~ scale_y_continuous(breaks = c(seq(0, 2.4, 0.60)), limits = c(0, 2.4), name = '', expand = c(0, 0)))) +
  geom_signif(y_position = c(2.10, 1.15), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p = 1', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(1.50, 1.2), xmin = c(1.8, 1.8), xmax = c(2.2, 2.2), annotation = c('p = 0.753', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2) +
  geom_signif(y_position = c(0.60, 1.15), xmin = c(2.8, 2.8), xmax = c(3.2, 3.2), annotation = c('p = 1', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2)


fig1c


fig1 <- cowplot::plot_grid(fig1a, fig1b, fig1c, nrow = 1)
fig1

#ggsave('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Figures/manuscript/NatureMedicine/Revision/Supplementary_Figure2_2025-05-09.pdf', fig1, width = 160, height = 40, units = 'mm', device = cairo_pdf)

###############################################################################
###############################################################################
### now for AMR only
################################################################################
################################################################################
################################################################################
################################################################################
## LOEUF
fig1a.amr <- ggplot(fig.df[fig.df$Pop %in% c('AMR') & fig.df$Variant_Type %in% c('PTV.0', 'PTV.1'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('PTV.0' = 'PTV 1\u02e2\u1d57-3\u02b3\u1d48\nLOEUF deciles',
                                                                 'PTV.1' = 'PTV 4\u1d57\u02b0-10\u1d57\u02b0\nLOEUF deciles')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = '') + ylab('Rate') +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 5, family = 'sans'), strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0, 0.01, -1, 0.01, 'cm'), panel.spacing = unit(0.1, 'lines'), legend.position = c(0.92, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0), size = 5, family = 'sans'), axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 12, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'PTV.0' ~ scale_y_continuous(breaks = c(seq(0, 0.3, 0.075)), limits = c(0, 0.3), expand = c(0, 0)),
    Variant_Type == 'PTV.1' ~ scale_y_continuous(breaks = c(seq(0, 0.3, 0.075)), limits = c(0, 0.3), expand = c(0, 0)))) +
  geom_signif(y_position = c(0.288, 0.27), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p < 5e-4', 'p = 0.753'), tip_length = 0, textsize = 0.9, size = 0.2)


fig1a.amr

## missense
fig1b.amr <- ggplot(fig.df[fig.df$Pop %in% c('AMR') & fig.df$Variant_Type %in% c('MisB', 'MisA'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = FALSE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('MisB' = 'missense\n MPC \u2265 2',
                                                                 'MisA' = 'missense\n 1 \u2264 MPC < 2')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = '') +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 5, family = 'sans'), 
        strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0, 0.01, -1, 0.01, 'cm'), panel.spacing = unit(0.1, 'lines'), 
        legend.position = c(0.92, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0), size = 5, family = 'sans'), 
        axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 12, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'MisB' ~ scale_y_continuous(breaks = c(seq(0, 0.64, 0.16)), limits = c(0, 0.64), name = '', expand = c(0, 0)),
    Variant_Type == 'MisA' ~ scale_y_continuous(breaks = c(seq(0, 0.64, 0.16)), limits = c(0, 0.64), name = '', expand = c(0, 0)))) +
  geom_signif(y_position = c(0.59, 0.235), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p = 0.077', 'p = 0.61'), tip_length = 0, textsize = 0.9, size = 0.2)


fig1b.amr

## synonymous
fig1c.amr <- ggplot(fig.df[fig.df$Pop %in% c('AMR') & fig.df$Variant_Type %in% c('MisX', 'SYN'),], aes(fill = Affected_Status, y = Rate_normalized, x = Pop)) +
  geom_bar(stat = 'identity', color = 'black', show.legend = TRUE, size = 0.2, width = 0.8, position = position_dodge()) +
  facet_wrap(~Variant_Type, labeller = labeller(Variant_Type = c('MisX' = 'missense\n MPC < 1',
                                                                 'SYN' = 'synonymous\n')), ncol = 2) + 
  theme_classic() +
  # labs(tag = expression(bold('a'))) +
  geom_errorbar(aes(x = Pop, ymin = CI.low.norm.adj, ymax = CI.high.norm.adj), position = position_dodge(0.8), width = 0.1, size = 0.2) +
  scale_fill_manual(values = c('#0072B2', '#56B4E9')) +
  scale_x_discrete(name = '') +
  theme(legend.title = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        legend.key.size = unit(0.2, 'cm'), legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0.1, 'cm'), 
        axis.text = element_text(size = 5, family = 'sans', color = 'black'), title = element_text(size = 5, family = 'sans'), strip.text.x = element_text(size = 5, family = 'sans', hjust = 0.5), strip.background = element_blank(),
        plot.margin = ggplot2::margin(0, 0.01, -1, 0.01, 'cm'), legend.position = c(0.80, 0.9), legend.background = element_blank(), legend.text = element_text(size = 5, family = 'sans'),
        axis.title.y = element_text(margin = ggplot2::margin(t = 0, r = 2, b = 0, l = 0), size = 5, family = 'sans'), axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 12, l = 0), size = 5, family = 'sans'), 
        axis.ticks.length = unit(.15, 'cm'), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.y = element_line(size = 0.2), axis.line.x = element_blank(), axis.ticks = element_line(size = 0.2)) +
  ggh4x::facetted_pos_scales(y = list(
    Variant_Type == 'MisX' ~ scale_y_continuous(breaks = c(seq(0, 2.36, 0.59)), limits = c(0, 2.36), name = '', expand = c(0, 0)),
    Variant_Type == 'SYN' ~ scale_y_continuous(breaks = c(seq(0, 2.36, 0.59)), limits = c(0, 2.36), name = '', expand = c(0, 0)))) +
  geom_signif(y_position = c(1.7, 1.2), xmin = c(0.8, 0.8), xmax = c(1.2, 1.2), annotation = c('p = 0.753', 'p = 1'), tip_length = 0, textsize = 0.9, size = 0.2)


fig1c.amr


fig1.amr <- cowplot::plot_grid(fig1a.amr, fig1b.amr, fig1c.amr, nrow = 1)
fig1.amr

# ggsave('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Figures/manuscript/NatureMedicine/Revision/Figure2_2025-05-09.pdf', fig1.amr, width = 120, height = 40, units = 'mm', device = cairo_pdf)

################################################################################
################################################################################
################################################################################
