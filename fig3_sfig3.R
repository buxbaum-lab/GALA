# Marina Natividad Avila
# 2024-12-09
# Fig 3 for Nature Medicine after Catalina's comments

# Let's fit linear models for Constraints Metrics projects
library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(cowplot)
library(data.table)


rm(list=ls())
`%notin%` <- Negate(`%in%`)

##################
##################
## data processing for ALL ancestries, at population maximum
afr <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/final/constraint_metrics_afr_downsample_8128.txt')
sas <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/final/constraint_metrics_sas_downsample_15308.txt')
eas <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/final/constraint_metrics_eas_downsample_9197.txt')
amr <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/final/constraint_metrics_amr_downsample_17296.txt')
nfe <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/final/constraint_metrics_nfe_downsample_56885.txt')

process_df <- function(df, pop){
  df$key <- paste(df$gene, df$transcript, sep = '::')
  df$canonical <- NULL; df$pop <- NULL; df$downsampling <- NULL; df$n_sites <- NULL; 
  df$gene <- NULL; df$transcript <- NULL; df$caf <- NULL; df$obs_lof <- NULL
  names(df)[names(df) == 'exp_syn'] <- paste('exp_syn', pop, sep = '_')
  names(df)[names(df) == 'obs_syn'] <- paste('obs_syn', pop, sep = '_')
  names(df)[names(df) == 'exp_mis'] <- paste('exp_mis', pop, sep = '_')
  names(df)[names(df) == 'obs_mis'] <- paste('obs_mis', pop, sep = '_')
  names(df)[names(df) == 'exp_lof'] <- paste('exp_lof', pop, sep = '_')
  names(df)[names(df) == 'obs_lof2'] <- paste('obs_lof', pop, sep = '_')
  #df$pop <- pop
  return(df)
}

afr <- process_df(afr, 'afr')
sas <- process_df(sas, 'sas')
eas <- process_df(eas, 'eas')
amr <- process_df(amr, 'amr')
nfe <- process_df(nfe, 'nfe')

max.pops <- merge(nfe, amr, by = 'key')
max.pops <- merge(max.pops, eas, by = 'key')
max.pops <- merge(max.pops, sas, by = 'key')
max.pops <- merge(max.pops, afr, by = 'key')

metrics <- fread("/sc/arion/projects/buxbaj01a/GALA/DATA/final_data/gnomAD/gnomad.v2.1.1.lof_metrics.by_gene.txt")
metrics$key <- paste(metrics$gene, metrics$transcript, sep = '::')
metrics <- metrics[, .(key, gene_id, chromosome, gene_length, gene_type, exp_syn, obs_syn, exp_mis, obs_mis, exp_lof, obs_lof, oe_lof_upper_bin)]

process_metrics <- function(df, pop){
  names(df)[names(df) == 'exp_syn'] <- paste('exp_syn', pop, sep = '_')
  names(df)[names(df) == 'obs_syn'] <- paste('obs_syn', pop, sep = '_')
  names(df)[names(df) == 'exp_mis'] <- paste('exp_mis', pop, sep = '_')
  names(df)[names(df) == 'obs_mis'] <- paste('obs_mis', pop, sep = '_')
  names(df)[names(df) == 'exp_lof'] <- paste('exp_lof', pop, sep = '_')
  names(df)[names(df) == 'obs_lof'] <- paste('obs_lof', pop, sep = '_')
  return(df)
}

metrics <- process_metrics(metrics, 'gnomAD')

df <- merge(metrics, max.pops, by = 'key')
next.df <- setDT(df)

# make counts df
deciles <- c(0:9)
df2 <- as.data.frame(deciles)
df2$length_genes <- NA # adding gene length based on gnomAD metrics$gene_length
df2$pLoF_AFR <- NA  
df2$pLoF_AMR <- NA  
df2$pLoF_EAS <- NA  
df2$pLoF_SAS <- NA  
df2$pLoF_NFE <- NA  

for (i in c(0:9)){
  tmp <- df[df$oe_lof_upper_bin == i]
  length_genes <- sum(tmp$gene_length, na.rm = TRUE)
  pLoF_afr <- sum(tmp$obs_lof_afr, na.rm = TRUE)
  pLoF_amr <- sum(tmp$obs_lof_amr, na.rm = TRUE)
  pLoF_nfe <- sum(tmp$obs_lof_nfe, na.rm = TRUE)
  pLoF_sas <- sum(tmp$obs_lof_sas, na.rm = TRUE)
  pLoF_eas <- sum(tmp$obs_lof_eas, na.rm = TRUE)
  df2$length_genes[df2$deciles == i] <- length_genes
  df2$pLoF_AFR[df2$deciles == i] <- pLoF_afr
  df2$pLoF_AMR[df2$deciles == i] <- pLoF_amr
  df2$pLoF_NFE[df2$deciles == i] <- pLoF_nfe
  df2$pLoF_SAS[df2$deciles == i] <- pLoF_sas
  df2$pLoF_EAS[df2$deciles == i] <- pLoF_eas
}

#####
# check gene lengths
# sum(metrics$gene_length) # 1,312,122,580
# sum(df2$length_genes) # 1,308,111,054
# dim(metrics[is.na(metrics$oe_lof_upper_bin)]) # 507  12... there's 507 genes with no LOEUF score (reasons explained in paper)
# 
# sum(metrics$gene_length[!is.na(metrics$oe_lof_upper_bin)]) == sum(df2$length_genes) # TRUE
# sum(metrics$gene_length[metrics$oe_lof_upper_bin == 0], na.rm = T) == sum(df2$length_genes[df2$deciles == 0])
# sum(metrics$gene_length[metrics$oe_lof_upper_bin == 5], na.rm = T) == sum(df2$length_genes[df2$deciles == 5])
# sum(metrics$gene_length[metrics$oe_lof_upper_bin == 9], na.rm = T) == sum(df2$length_genes[df2$deciles == 9])
#####

# total length over all genes
exome.size <- sum(df2$length_genes) # 1308111054

###
scaled_df <- setDT(df2)


#Probably should normalize the numerator and denominator of the ratio to sample size. 
n.afr <- 8128; n.amr <- 17296; n.eas <- 9197; n.sas <- 15308; n.nfe <- 56885

# Note. The comment below was from the email where Joseph asked to scale.
# after talking to Xuran, this is how it should be scaled
scaled_df$AFR <- scaled_df$pLoF_AFR #(scaled_df$pLoF_AFR/n.afr)/(scaled_df$pLoF_NFE/n.nfe)
scaled_df$AMR <- scaled_df$pLoF_AMR #(scaled_df$pLoF_AMR/n.amr)/(scaled_df$pLoF_NFE/n.nfe)
scaled_df$EAS <- scaled_df$pLoF_EAS #(scaled_df$pLoF_EAS/n.eas)/(scaled_df$pLoF_NFE/n.nfe)
scaled_df$SAS <- scaled_df$pLoF_SAS #(scaled_df$pLoF_SAS/n.sas)/(scaled_df$pLoF_NFE/n.nfe)
scaled_df$NFE <- scaled_df$pLoF_NFE #(scaled_df$pLoF_SAS/n.sas)/(scaled_df$pLoF_NFE/n.nfe)


#melt data frame into long format
tst <- reshape2::melt(scaled_df ,  id.vars = 'deciles', variable.name = 'pop')
tst <- tst[tst$pop %in% c('AMR', 'AFR', 'EAS', 'SAS','NFE'),]

tst$afr <- 0
tst$afr[tst$pop == 'AFR'] <- 1
tst$amr <- 0
tst$amr[tst$pop == 'AMR'] <- 1
tst$sas <- 0
tst$sas[tst$pop == 'SAS'] <- 1
tst$eas <- 0
tst$eas[tst$pop == 'EAS'] <- 1
tst$nfe <- 0
tst$nfe[tst$pop == 'NFE'] <- 1

tst$pop <- NULL
tst$pop[tst$afr == 1] <- 'AFR'
tst$pop[tst$amr == 1] <- 'AMR'
tst$pop[tst$sas == 1] <- 'SAS'
tst$pop[tst$eas == 1] <- 'EAS'
tst$pop[tst$nfe == 1] <- 'NFE'


# For the legend
tst$Population <- NA
tst$Population[tst$pop == 'AFR'] <- 'AFR: 8,128'
tst$Population[tst$pop == 'AMR'] <- 'AMR: 17,296'
tst$Population[tst$pop == 'EAS'] <- 'EAS: 9,197'
tst$Population[tst$pop == 'SAS'] <- 'SAS: 15,308'
tst$Population[tst$pop == 'NFE'] <- 'NFE: 56,885'

# sample size
tst$n.samples <- NA
tst$n.samples[tst$pop == 'AFR'] <- n.afr
tst$n.samples[tst$pop == 'AMR'] <- n.amr
tst$n.samples[tst$pop == 'EAS'] <- n.eas
tst$n.samples[tst$pop == 'SAS'] <- n.sas
tst$n.samples[tst$pop == 'NFE'] <- n.nfe

# gene length 
tst$gene.length <- NA
tst$gene.length <- df2$length_genes[match(tst$deciles, df2$deciles)]

# scaled gene length
tst$scaled.gene.length <- tst$gene.length/exome.size

names(tst)[names(tst) == 'deciles'] <- 'constraint'

rownames(tst) <- paste(tst$pop, '::', 'bin_', tst$constraint, sep = '')

# write.table(tst, '/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Constraint/table_for_lm_2023-10-10.txt', col.names = T, row.names = T, quote = F, sep = '\t')

# scale over population size
tst$value2 <- tst$value/tst$n.samples
# scale over gene length
tst$value3 <- tst$value2/tst$scaled.gene.length

###########################
###########################

# plot!
b <- tst %>% ggplot + 
  aes(x = constraint, y = value3) + geom_point(aes(fill = factor(pop), shape = factor(pop)), colour = 'black', size = 4) + 
  theme_classic() +
  labs(x = 'gnomAD LOEUF',
       y = expression(paste('(\u03A3 observed ', PTV['pop'], ') / (', N['pop'], 'x gene length)')),
       title = '') +
       #subtitle = expression(italic('n.genes = 19,197')),
  scale_y_continuous(breaks = seq(0, 10, 1.0), limits=c(0, 9.1)) +
  scale_size(guide = 'none') +
  theme(legend.title = element_blank(), axis.text = element_text(size = 5, color = 'black'),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), element_blank(),
        axis.ticks.length = unit(.15, 'cm'), plot.margin = ggplot2::margin(0.025, 0.025, 0.025, 0.025, 'cm'), 
        axis.title.y = element_text(margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1), size = 5), 
        axis.title = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 5),
        axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2),
        legend.position = c(0.15, 0.79), legend.direction = 'vertical', legend.text = element_text(size = 5)) + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_shape_manual(values=c(22, 21, 24, 23, 25)) +
  scale_fill_manual(values = c('#F8766D', '#A2A500', '#00BF7D', '#00AFF6', '#E76BF3')) +
  scale_x_continuous(breaks = 0:9, #expand = c(0.05, 0),
                     labels = c(expression(1), expression(2), expression(3), expression(4), expression(5), 
                                expression(6), expression(7), expression(8), expression(9), expression(10)))


b

##############
##############
# only AMR and NFE
# plot!
amr.nfe <-  tst[tst$pop %in% c('AMR', 'NFE'),] %>% ggplot + 
  aes(x = constraint, y = value3) + geom_point(aes(fill = factor(pop), shape = factor(pop)), colour = 'black', size = 4) + 
  theme_classic() +
  labs(x = 'gnomAD LOEUF',
       y = expression(paste('(\u03A3 observed ', PTV['pop'], ') / (', N['pop'], 'x gene length)')),
       title = '') +
  scale_y_continuous(breaks = seq(0, 6, 1), limits=c(0, 6)) +
  scale_size(guide = 'none') +
  theme(legend.title = element_blank(), axis.text = element_text(size = 5, color = 'black'),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), element_blank(),
        axis.ticks.length = unit(.25, 'cm'), plot.margin = ggplot2::margin(0.025, 0.025, 0.025, 0.025, 'cm'), 
        axis.title.y = element_text(margin = ggplot2::margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1), size = 5), 
        axis.title = element_text(margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), size = 5),
        axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2),
        legend.position = c(0.9, 0.15), legend.direction = 'vertical', legend.text = element_text(size = 5)) + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_shape_manual(values=c(21, 23)) +
  scale_fill_manual(values = c('#A3A500', '#00B0F6')) +
  scale_x_continuous(breaks = 0:9, #expand = c(0.05, 0),
                     labels = c(expression(1), expression(2), expression(3), expression(4), expression(5), 
                                expression(6), expression(7), expression(8), expression(9), expression(10)))


amr.nfe

# ggsave('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/Figures/manuscript/NatureMedicine/new/Figure3_2024-11-27.pdf', amr.nfe, height = 80, width = 80, units = 'mm', device = cairo_pdf)
