# This script was adapted from code made available by the Talkowski Lab:
# https://github.com/talkowski-lab/TADA_2022
# Marina Natividad Avila
# 2024-04-17

library(GenomicRanges); library(stringr); library(openxlsx)
rm(list=ls())
`%notin%` <- Negate(`%in%`)

### Change to working directory where the folder was  
source('/sc/arion/projects/buxbaj01a/software/TADA_fu_prerelease/functions.R')
### File name for supplementary table
supp_tab <- '/sc/arion/projects/buxbaj01a/GALA/DATA/final_data/Fu/41588_2022_1104_MOESM3_ESM.xlsx'

########################################################################
### List of ACMG genes that were not released in certain cohorts
########################################################################
acmg <- c('ACTA2', 'ACTC1', 'APC', 'APOB', 'ATP7B', 'BMPR1A', 'BRCA1', 'BRCA2', 'CACNA1S', 'COL3A1', 'DSC2', 'DSG2', 'DSP', 'FBN1', 'GLA', 'KCNH2', 'KCNQ1', 'LDLR', 'LMNA', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'MYBPC3', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NF2', 'OTC', 'PCSK9', 'PKP2', 'PMS2', 'PRKAG2', 'RB1', 'RET', 'RYR1', 'RYR2', 'SCN5A', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD3', 'SMAD4', 'STK11', 'TGFBR1', 'TGFBR2', 'TMEM43', 'TNNI3', 'TNNT2', 'TP53', 'TPM1', 'VHL', 'WT1')

### Information on genes: constraint, priors
info <- openxlsx::read.xlsx(supp_tab, sheet = 'Supplementary Table 8')
info <- info[1:18128,]

final <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/manuscript/2024/AMR_OTH_pedigree_2024_04_03.txt')
final$isAMR <- 0
final$isAMR[final$Pop == 'AMR'] <- 1
with(final, table(Affected_Status, Pop))
#                 Pop
# Affected_Status   AFR   AMR   EAS   FIN   MID   NFE   SAS
#               1   184  1459   112     7   121  4079   246
#               2   477  4450   419   117   446 10998   573

with(final, table(isAMR, Affected_Status))
#      Affected_Status
# isAMR     1     2
#     0  4749 13030
#     1  1459  4450

##############################
##############################
# now cnv numbers
cnv_ped <- final[final$CNV_called == TRUE & final$CNV_called_father == TRUE & final$CNV_called_mother == TRUE]
with(cnv_ped, table(isAMR, Affected_Status))
#   Affected_Status
# isAMR     1     2
#     0  4186 11321
#     1   705  2152
##############################
###  Setting sample sizes  ###
##############################
## ASC fu
n_prob_asc <- 963; n_sib_asc <- 325
## GALA
n_prob_asn <- 440; n_sib_asn <- 26

## SPARK Fu
n_prob_spk <- 1041; n_sib_spk <- 417
## SPARK iWESv1 and iWESv2
n_prob_spp <- 2006; n_sib_spp <- 691
## Case control
n1_case = 267; n1_control = 801 # BioMe controls

## CNV
n_prob_cnv <- 2152; n_sib_cnv <- 705

# n_prob_inh_cnv <- 330; n_sib_inh_cnv <- 27
n_case_cnv <- 250; n_control_cnv <- 752

### DD numbers
n_prob_dd <- 31058; n_sin_dd <- 0

###################################
###  Loading gene count tables  ###
###################################
### SNV, denovo + inherited, ASC+GALA (sequenced at Broad Institute)
dg_asc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/merged-tables/amr_asc_fu_2023-11-06.txt')
dg_gal <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/2024/GALA_denovo_inherited_2024-04-04.txt')

### SNV, denovo + inherited, SPARK (sequenced at Regeneron)
dg_sp0 <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/merged-tables/amr_spark_fu_2023-11-06.txt')
dg_sp1 <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/2024/amr_SPARK_iWESv2_denovo_inherited_2024-04-04.txt')

# case control counts
dg_cc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/case-control_2024-03-20.txt')

dg_ddd <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/dg_ddd_2024-02-08.txt')

## parameters ##
S <- sort(info$gene_id,index.return=TRUE)
  info <- info[S$ix,]
S <- sort(dg_asc$Gene_ID,index.return=TRUE)
  dg_asc <- dg_asc[S$ix,]
S <- sort(dg_gal$Gene_ID,index.return=TRUE)
  dg_gal <- dg_gal[S$ix,]

S <- sort(dg_sp0$Gene_ID,index.return=TRUE)
  dg_sp0 <- dg_sp0[S$ix,]
S <- sort(dg_sp1$Gene_ID,index.return=TRUE)
  dg_sp1 <- dg_sp1[S$ix,]

S <- sort(dg_cc$ENSG,index.return=TRUE)
  dg_cc <- dg_cc[S$ix,]

S <- sort(dg_ddd$Gene_ID,index.return=TRUE)
  dg_ddd <- dg_ddd[S$ix,]

##################################
###  calibrate mutation rates  ###
##################################
# BI-sequenced
cal_asc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/calibration-tables/calibration-asc_2024-03-29.txt')

  dg_asc$mut.ptv <- cal_asc$mut.ptv[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]
  dg_asc$mut.misb <- cal_asc$mut.misb[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]
  dg_asc$mut.misa <- cal_asc$mut.misa[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]

# Regeneron-sequenced
cal_spk <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/calibration-tables/calibration-spk_2024-03-29.txt')

  dg_sp0$mut.ptv <- cal_spk$mut.ptv[match(dg_sp0$Gene_ID, cal_spk$Gene_ID)]
  dg_sp0$mut.misb <- cal_spk$mut.misb[match(dg_sp0$Gene_ID, cal_spk$Gene_ID)]
  dg_sp0$mut.misa <- cal_spk$mut.misa[match(dg_sp0$Gene_ID, cal_spk$Gene_ID)]

#################################
### Loading CNV count tables  ###
#################################

### CNV, denovo, ASC+SSC+SPARK
cnv_dn <- fread('/sc/arion/projects/buxbaj01a/GALA/DATA/edited_data/TADA_August2023/AMR_ASC_CNV_counts_2023-09-08.txt')
cnv_dn <- cnv_dn[-which(cnv_dn$chr=="chr21" & (cnv_dn$end-cnv_dn$start)>30000000),] ### remove trisomy 21 94 --> 87
cnv_dn <- cnv_dn[cnv_dn$non_diploid_freq < 0.01,] # 87 --> 87
cnv_dn$notes <- NULL
cnv_dn$non_diploid_freq <- NULL
cnv_dn$non_diploid_count <- NULL
cnv_dn$isAMR <- NULL
  cnv_ped <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pedigree/November2023/GALA_AMR_ped_2023-11-29.txt')
  cnv_dn <- cnv_dn[cnv_dn$sample %in% ped$Sample]

cnv_dn2 <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/CNV/dn_prospect_annotated.txt')
cnv_dn2 <- cnv_dn2[-which(cnv_dn2$chr==21 & (cnv_dn2$end-cnv_dn2$start)>30000000),] # 140 --> 139
cnv_dn2$chr <- paste('chr', cnv_dn2$chr, sep = '')

names(cnv_dn2)[names(cnv_dn2) == 'svtype'] <- 'call'
names(cnv_dn2)[names(cnv_dn2) == 'vaf'] <- 'site_freq'
names(cnv_dn2)[names(cnv_dn2) == 'batch'] <- 'cluster'
names(cnv_dn2)[names(cnv_dn2) == 'variant_name'] <- 'site_name'
names(cnv_dn2)[names(cnv_dn2) == 'vac'] <- 'site_count'
# names(cnv_dn2)[names(cnv_dn2) == 'genes_strict_overlap_totalExons'] <- 'num_exons'
names(cnv_dn2)[names(cnv_dn2) == 'genes_any_overlap_totalExons'] <- 'num_exons'
names(cnv_dn2)[names(cnv_dn2) == 'known_CNV_gd_id'] <- 'gd_loci'
names(cnv_dn2)[names(cnv_dn2) == 'known_CNV_nahr'] <- 'nahr'
#names(cnv_dn2)[names(cnv_dn2) == 'genes_strict_overlap'] <- 'genes'
names(cnv_dn2)[names(cnv_dn2) == 'genes_any_overlap'] <- 'genes'

cnv_dn2$Dataset <- 'GALA'
  to_keep <- names(cnv_dn)

tmp <- cnv_dn2[, c('chr', 'start', 'end', 'QS', 'CN', 'call', 'sample', 'cluster', 'site_name', 'site_count',
                   'site_freq', 'num_exons', 'Affected_Status', 'gd_loci', 'nahr', 'genes', 'Dataset')]

tmp$genes[tmp$genes == 'None'] <- NA

cnv_dn <- rbind(cnv_dn, tmp)
cnv_dn$nahr[is.na(cnv_dn$nahr)] <- FALSE
table(cnv_dn$nahr)
# FALSE  TRUE 
#   188    38 

### CNV case-control, ASC
cnv_cc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/cases-BioMe/CNV/GALA_BioMe_GATK_gCNV_250cases_and_752_controls_685CNVs_input_GD_fixed.txt')
# cnv_cc <- cnv_cc[-which(cnv_cc$chr == 21 & (cnv_cc$end-cnv_cc$start)>30000000)] ### no trisomy 21 present
cnv_cc$chr <- paste('chr', cnv_cc$chr, sep = '')
cnv_cc$Affected_Status <- NA
cnv_cc$Affected_Status[cnv_cc$data == 'case'] <- 2
cnv_cc$Affected_Status[cnv_cc$data == 'control'] <- 1
dim(cnv_cc[is.na(cnv_cc$Affected_Status)])

cnv_cc$genes[cnv_cc$genes == 'None'] <- NA
cnv_cc$nahr[is.na(cnv_cc$nahr)] <- FALSE
cnv_cc$GD_loci <- FALSE
cnv_cc$GD_loci[!is.na(cnv_cc$c22.gd_id)] <- TRUE
with(cnv_cc, table(nahr, GD_loci))
#       GD_loci
# nahr    FALSE TRUE
#   FALSE   676    0
#   TRUE      0    9

###################################
###  Bayes Factor calculations  ###
###################################
beta.dn <- 0.2

###################
###  SNV/indel  ###
###################

###################
### de novo
### DN PTV
BF_dn_ptv_asc <- BF_DN_SNV(count_case=dg_asc$dn.ptv+dg_gal$dn.ptv, count_con=dg_asc$dn.ptv.sib+dg_gal$dn.ptv.sib, n_case=n_prob_asc+n_prob_asn, n_con=n_sib_asc+n_sib_asn, mut=dg_asc$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv_spk <- BF_DN_SNV(count_case=dg_sp0$dn.ptv+dg_sp1$dn.ptv, count_con=dg_sp0$dn.ptv.sib+dg_sp1$dn.ptv.sib, n_case=n_prob_spk+n_prob_spp, n_con=n_sib_spk+n_sib_spp, mut=dg_sp0$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

### DN misB
BF_dn_misB_asc <- BF_DN_SNV(count_case=dg_asc$dn.misb+dg_gal$dn.misb, count_con=dg_asc$dn.misb.sib+dg_gal$dn.misb.sib, n_case=n_prob_asc+n_prob_asn, n_con=n_sib_asc+n_sib_asn, mut=dg_asc$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB_spk <- BF_DN_SNV(count_case=dg_sp0$dn.misb+dg_sp1$dn.misb, count_con=dg_sp0$dn.misb.sib+dg_sp1$dn.misb.sib, n_case=n_prob_spk+n_prob_spp, n_con=n_sib_spk+n_sib_spp, mut=dg_sp0$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

### DN misA
BF_dn_misA_asc <- BF_DN_SNV(count_case=dg_asc$dn.misa+dg_gal$dn.misa, count_con=dg_asc$dn.misa.sib+dg_gal$dn.misa.sib, n_case=n_prob_asc+n_prob_asn, n_con=n_sib_asc+n_sib_asn, mut=dg_asc$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA_spk <- BF_DN_SNV(count_case=dg_sp0$dn.misa+dg_sp1$dn.misa, count_con=dg_sp0$dn.misa.sib+dg_sp1$dn.misa.sib, n_case=n_prob_spk+n_prob_spp, n_con=n_sib_spk+n_sib_spp, mut=dg_sp0$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk

###################
### inherited
### IN PTV
BF_in_ptv <- BF_CC_SNV(count_case=dg_asc$t.ptv+dg_gal$t.ptv+dg_sp0$t.ptv+dg_sp1$t.ptv, count_con=dg_asc$u.ptv+dg_gal$u.ptv+dg_sp0$u.ptv+dg_sp1$u.ptv, n_case=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, n_con=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, mut=dg_asc$mut.ptv, prior=info$prior.in.ptv)
### IN missenseB
BF_in_misB <- BF_CC_SNV(count_case=dg_asc$t.misb+dg_gal$t.misb+dg_sp0$t.misb+dg_sp1$t.misb, count_con=dg_asc$u.misb+dg_gal$u.misb+dg_sp0$u.misb+dg_sp1$u.misb,  n_case=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, n_con=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, mut=dg_asc$mut.misb, prior=info$prior.in.misb)
### IN missenseA
BF_in_misA <- BF_CC_SNV(count_case=dg_asc$t.misa+dg_gal$t.misa+dg_sp0$t.misa+dg_sp1$t.misa, count_con=dg_asc$u.misa+dg_gal$u.misa+dg_sp0$u.misa+dg_sp1$u.misa,  n_case=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, n_con=n_prob_asc+n_prob_asn+n_prob_spk+n_prob_spp, mut=dg_asc$mut.misa, prior=info$prior.in.misa)


###################
### case-control
### CC PTV
  BF_cc_ptv <- BF_CC_SNV(count_case=dg_cc$case.ptv, count_con=dg_cc$ctl.ptv, n_case=n1_case, n_con=n1_control, mut=dg_asc$mut.ptv, prior=info$prior.cc.ptv)
### CC misB
  BF_cc_misB <- BF_CC_SNV(count_case=dg_cc$case.misb, count_con=dg_cc$ctl.misb, n_case=n1_case, n_con=n1_control, mut=dg_asc$mut.misb, prior=info$prior.cc.misb)
### CC misA
  BF_cc_misA <- BF_CC_SNV(count_case=dg_cc$case.misa, count_con=dg_cc$ctl.misa, n_case=n1_case, n_con=n1_control, mut=dg_asc$mut.misa, prior=info$prior.cc.misa)

############
### CNV  ###
############
# size range, in number of constrained genes, to consider
cnv_size_range <- 1:8 

###################
### de novo
cnv <- processCNV(cnv_dn, loeuf_threshold=0.6,  info=info)
# cnv_use <- cnv[which((is.na(cnv$gd_loci) | cnv$nahr==FALSE) & cnv$non_diploid_freq<.01 & cnv$num_genes %in% cnv_size_range & !is.na(cnv$Affected_Status)),]

table(cnv$num_genes)
#  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 20 
# 98 59 14  5  6  3  3  6 11  2  2  3  1  2  7  2  1  1 

cnv_use <- cnv[which((cnv$gd_loci == FALSE | cnv$nahr==FALSE) & cnv$num_genes %in% cnv_size_range & !is.na(cnv$Affected_Status)),] # 79
  cnv_use$key <- paste(cnv_use$site_name, cnv_use$sample, sep = '::')
  cnv_use <- cnv_use[cnv_use$key %notin% c('suffix_27472::COL-065C', 'suffix_20698::cc1538.201', 'suffix_27472::COL-116C')] # 79 --> 77
  cnv_use$key <- NULL
  
del_dup_adj <- table(cnv[,c("Affected_Status", "call")])
del_dup_adj <- (del_dup_adj[2,1]/del_dup_adj[1,1])/(del_dup_adj[2,2]/del_dup_adj[1,2])

del_use <- cnv_use[cnv_use$call=="DEL",]
mut.pred.del <- sapply(cnv_size_range, function(x){ max(1, sum(del_use$num_genes==x &del_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
del_use$mut <- mut.pred.del[del_use$num_genes]

BF_dn_del_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv)
BF_dn_del_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv)
BF_dn_del <- pmax(1, BF_dn_del_prob/BF_dn_del_sib)

dup_use <- cnv_use[cnv_use$call=="DUP",]
mut.pred.dup <- sapply(cnv_size_range, function(x){ max(1, sum(dup_use$num_genes==x &dup_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
dup_use$mut <- mut.pred.dup[dup_use$num_genes]

BF_dn_dup_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv, del_dup_adj=del_dup_adj)
BF_dn_dup_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv, del_dup_adj=del_dup_adj)
BF_dn_dup <- pmax(1, BF_dn_dup_prob/BF_dn_dup_sib)

########
## CC
cnv_size_range <- 1:8
 
# cnv <- processCNV(cnv_cc[is.na(cnv_cc$c22.gd_id) | cnv_cc$nahr == FALSE], loeuf_threshold = 0.6, info = info)
cnv <- processCNV(cnv_cc[which((is.na(cnv_cc$c22.gd_id) & is.na(cnv_cc$c22.gd_id.1)) | cnv_cc$nahr==FALSE)], loeuf_threshold=0.6, info=info)
cnv_use <- cnv[cnv$num_genes %in% cnv_size_range,] # 648 --> 158
cnv_use$key <- paste(cnv_use$s, cnv_use$chr, sep = '::')
cnv_use <- cnv_use[cnv_use$key %notin% c('C1468_AU214A_v1_Exome_GCP::chr11', 'RP-1944_AFI2821446A_v1_Exome_GCP::chr5', 'RP-1944_AFI1960729A_v1_Exome_GCP::chr17')] # 160

BF_cc_del <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DEL",], n_case=n_case_cnv, n_con=n_control_cnv, info=info, prior=info$prior.cc.ptv, del_dup_adj=1, nu=5000)
BF_cc_dup <- BF_CC_CNV(cnv_size_range, cnv_use=cnv_use[cnv_use$call=="DUP",], n_case=n_case_cnv, n_con=n_control_cnv, info=info, prior=info$prior.cc.ptv, del_dup_adj=del_dup_adj, nu=5000)

################################################
###  Accounting for genes with ACMG masking  ###
################################################
acmg_ind <- which(info$gene %in% acmg)

BF_dn_ptv_spk_acmg <- BF_DN_SNV(count_case=dg_sp0$dn.ptv+dg_sp1$dn.ptv, count_con=dg_sp0$dn.ptv.sib+dg_sp1$dn.ptv.sib, n_case=n_prob_spp+n_prob_spk, n_con=0, mut=dg_sp0$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv_spk[acmg_ind] <- BF_dn_ptv_spk_acmg[acmg_ind]
BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

BF_dn_misB_spk_acmg <- BF_DN_SNV(count_case=dg_sp0$dn.misb+dg_sp1$dn.misb, count_con=dg_sp0$dn.misb.sib+dg_sp1$dn.misb.sib, n_case=n_prob_spp+n_prob_spk, n_con=0, mut=dg_sp0$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB_spk[acmg_ind] <- BF_dn_misB_spk_acmg[acmg_ind]
BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

BF_dn_misA_spk_acmg <- BF_DN_SNV(count_case=dg_sp0$dn.misa+dg_sp1$dn.misa, count_con=dg_sp0$dn.misa.sib+dg_sp1$dn.misa.sib, n_case=n_prob_spp+n_prob_spk, n_con=0, mut=dg_sp0$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA_spk[acmg_ind] <- BF_dn_misA_spk_acmg[acmg_ind]
BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk

########################################################################
### Integrating DDD - SNV/indel
########################################################################

### DN PTV
BF_dn_ptv_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.ptv, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)

### DN misB
BF_dn_misB_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misb, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)

### DN misA
BF_dn_misA_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misa, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)

##################################
###  TADA outcome calculation  ###
##################################

### Aggregate Bayes Factors
  # de novo only
BF_asd_dn <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk, 1),
                   pmax(BF_dn_misB_asc*BF_dn_misB_spk, 1),
                   pmax(BF_dn_misA_asc*BF_dn_misA_spk, 1),
                   pmax(BF_dn_del, 1), 
                   pmax(BF_dn_dup, 1))

  # de novo + inherited
BF_asd_inh <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv, 1),
                    pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB, 1),
                    pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA, 1),
                    pmax(BF_dn_del, 1), 
                    pmax(BF_dn_dup, 1))

  # de novo + inherited + case-control
BF_asd_cc <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv*BF_cc_ptv, 1),
                   pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB*BF_cc_misB, 1),
                   pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA*BF_cc_misA, 1),
                   pmax(BF_dn_del*BF_cc_del, 1), 
                   pmax(BF_dn_dup*BF_cc_dup, 1))


### Aggregate DD BF
BF_ddd <- cbind(pmax(1, BF_dn_ptv_ddd), pmax(BF_dn_misB_ddd, 1), pmax(BF_dn_misA_ddd, 1))

# Remove CNV evidence after aggregating 
### Aggregate NDD BF
BF_asd_ddd <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv*BF_cc_ptv*BF_dn_ptv_ddd, 1),
                    pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB*BF_cc_misB*BF_dn_misB_ddd, 1),
                    pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA*BF_cc_misA*BF_dn_misA_ddd, 1),
                    pmax(BF_dn_del*BF_cc_del, 1), 
                    pmax(BF_dn_dup*BF_cc_dup, 1))


 # Only use del/dup if there is snv/indel evidence (Tot BF > 2) 
# ASD evidence
BF_asd_cc[apply(BF_asd_cc[,1:3], 1, max)<5,4:5] <- 1
# ASD + DD evidence
BF_asd_ddd[apply(BF_asd_ddd[,1:3], 1, max)<5,4:5] <- 1

  
### FDR calculation
  # de novo only
qval_asd_dn <- Bayesian.FDR(apply(BF_asd_dn, 1, prod), pi0 = 1 - 0.05)
  names(qval_asd_dn) <- info$gene_id

  # de novo + inherited
qval_asd_inh <- Bayesian.FDR(apply(BF_asd_inh, 1, prod), pi0 = 1 - 0.05)
  names(qval_asd_inh) <- info$gene_id

  # de novo + inherited + case-control
qval_asd_cc <- Bayesian.FDR(apply(BF_asd_cc, 1, prod), pi0 = 1 - 0.05)
  names(qval_asd_cc) <- info$gene_id

  # DD
qval_ddd <- Bayesian.FDR(apply(BF_ddd, 1, prod), pi0 = 1 - 0.05)
  names(qval_ddd) <- info$gene_id
  
  # DD + GALA
qval_asd_ddd <- Bayesian.FDR(apply(BF_asd_ddd, 1, prod), pi0 = 1 - 0.05)
  names(qval_asd_ddd) <- info$gene_id
  
### TADA-ASD significant genes at various thresholds
sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_dn<x))
# 7 14 28 50     de novo only

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_inh<x))
# 7 15 33 58  de novo + inherited

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_cc<x))
# 8 16 35 61  de novo + inherited + case-control

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_ddd<x))
# 309 378 477 559 DD

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_ddd<x))
# 326 415 545 654


### back-transform FDR to p-values
  # de novo only
ec_threshold <- .05/length(qval_asd_dn)
pval_asd_dn <- q2p(qval_asd_dn)
sum(pval_asd_dn <= ec_threshold) # 10

# de novo + inherited
ec_threshold <- .05/length(qval_asd_inh)
pval_asd_inh <- q2p(qval_asd_inh)
sum(qval_asd_inh <= ec_threshold) # 6

# de novo + inherited + case-control
ec_threshold <- .05/length(qval_asd_cc)
pval_asd_cc <- q2p(qval_asd_cc)
sum(qval_asd_cc <= ec_threshold) # 4

# DD
ec_threshold <- .05/length(qval_ddd)
pval_ddd <- q2p(qval_ddd)
sum(qval_ddd <= ec_threshold) # 212

# DD + GALA
ec_threshold <- .05/length(qval_asd_ddd)
pval_asd_ddd <- q2p(qval_asd_ddd)
sum(qval_asd_ddd <= ec_threshold) # 227

#######
## do DD values make sense?
fu <- openxlsx::read.xlsx('/sc/arion/projects/buxbaj01a/GALA/DATA/final_data/Fu/41588_2022_1104_MOESM3_ESM.xlsx', sheet = 'Supplementary Table 11')
  fu <- fu[1:18128,]
  fu <- setDT(fu)

# numbers should be 309 378 477 559
dim(fu[fu$FDR_TADA_DD <= 0.001])  # 309
dim(fu[fu$FDR_TADA_DD <= 0.01])   # 378
dim(fu[fu$FDR_TADA_DD <= 0.05])   # 477
dim(fu[fu$FDR_TADA_DD <= 0.1])    # 559

# # de novo table
# BF_table <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk, 1),
#                   pmax(BF_dn_misB_asc*BF_dn_misB_spk, 1),
#                   pmax(BF_dn_misA_asc*BF_dn_misA_spk, 1),
#                   pmax(BF_dn_del, 1), pmax(BF_dn_dup, 1))

# de novo + inherited + case-control
BF_table <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk, 1),
                  pmax(BF_dn_misB_asc*BF_dn_misB_spk, 1),
                  pmax(BF_dn_misA_asc*BF_dn_misA_spk, 1),
                  pmax(BF_dn_del, 1), 
                  pmax(BF_dn_dup, 1),
                  pmax(BF_in_ptv, 1), 
                  pmax(BF_in_misB, 1), 
                  pmax(BF_in_misA, 1),
                  pmax(BF_cc_ptv, 1), 
                  pmax(BF_cc_misB, 1), 
                  pmax(BF_cc_misA, 1), 
                  pmax(BF_cc_del, 1), 
                  pmax(BF_cc_dup, 1),
                  pmax(BF_dn_ptv_ddd, 1),
                  pmax(BF_dn_misB_ddd, 1),
                  pmax(BF_dn_misA_ddd, 1))

df <- as.data.frame(BF_table)
df2 <- as.data.frame(BF_asd_cc)

names(df) <- c('BF.dn.ptv', 'BF.dn.misb', 'BF.dn.misa', 'BF.dn.del', 'BF.dn.dup', 
               'BF.in.ptv', 'BF.in.misb', 'BF.in.misa', 
               'BF.cc.ptv', 'BF.cc.misb', 'BF.cc.misa', 'BF.cc.del', 'BF.cc.dup',
               'BF.dn.ptv.ddd', 'BF.dn.misb.ddd', 'BF.dn.misa.ddd')

names(df2) <- c('BF.ptv', 'BF.misb', 'BF.misa', 'BF.del', 'BF.dup')


res <- data.table(gene = info$gene, gene_id = info$gene_id, pLI = info$pLI, LOEUF = info$LOEUF,
                  BF.dn.ptv = df$BF.dn.ptv, BF.dn.misb = df$BF.dn.misb, BF.dn.misa = df$BF.dn.misa, BF.dn.del = df$BF.dn.del, BF.dn.dup = df$BF.dn.dup,
                  BF.in.ptv = df$BF.in.ptv, BF.in.misB = df$BF.in.misb, BF.in.misA = df$BF.in.misa, 
                  BF.cc.ptv = df$BF.cc.ptv, BF.cc.misb = df$BF.cc.misb, BF.cc.misa = df$BF.cc.misa, BF.cc.del = df$BF.cc.del, BF.cc.dup = df$BF.cc.dup,
                  BF.dn.ptv.ddd = df$BF.dn.ptv.ddd, BF.dn.misb.ddd = df$BF.dn.misb.ddd, BF.dn.misa.ddd = df$BF.dn.misa.ddd,
                  qval_dn = qval_asd_dn, pval_dn = pval_asd_dn, qval_inh = qval_asd_inh, pval_inh = pval_asd_inh, qval_cc = qval_asd_cc, pval_cc = pval_asd_cc, 
                  qval_gala_dd = qval_asd_ddd, pval_gala_dd = pval_asd_ddd, qval_dd = qval_ddd, pval_dd = pval_ddd,
                  l10PTV = log10(df2$BF.ptv), l10misb = log10(df2$BF.misb), l10misa = log10(df2$BF.misa), l10DEL = log10(df2$BF.del), l10DUP = log10(df2$BF.dup))

# write.table(res, '/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/manuscript/2024/TADA_results/amr_bf_table_2024-04-17.txt', col.names = T, row.names = F, quote = F, sep = '\t')

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

## non-AMR Fu et al. 2022
rm(list=ls())
`%notin%` <- Negate(`%in%`)

### Change to working directory where the folder was cloned
source('/sc/arion/projects/buxbaj01a/software/TADA_fu_prerelease/functions.R')
### File name for supplementary table
supp_tab <- '/sc/arion/projects/buxbaj01a/GALA/DATA/final_data/Fu/41588_2022_1104_MOESM3_ESM.xlsx'

ped <- openxlsx::read.xlsx(supp_tab, sheet = 'Supplementary Table 4')

########################################################################
### List of ACMG genes that were not released in certain cohorts
########################################################################
acmg <- c('ACTA2', 'ACTC1', 'APC', 'APOB', 'ATP7B', 'BMPR1A', 'BRCA1', 'BRCA2', 'CACNA1S', 'COL3A1', 'DSC2', 'DSG2', 'DSP', 'FBN1', 'GLA', 'KCNH2', 'KCNQ1', 'LDLR', 'LMNA', 'MEN1', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'MYBPC3', 'MYH11', 'MYH7', 'MYL2', 'MYL3', 'NF2', 'OTC', 'PCSK9', 'PKP2', 'PMS2', 'PRKAG2', 'RB1', 'RET', 'RYR1', 'RYR2', 'SCN5A', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SMAD3', 'SMAD4', 'STK11', 'TGFBR1', 'TGFBR2', 'TMEM43', 'TNNI3', 'TNNT2', 'TP53', 'TPM1', 'VHL', 'WT1')

### Information on genes: constraint, priors
info <- openxlsx::read.xlsx(supp_tab, sheet = 'Supplementary Table 8')
info <- info[1:18128,]

##############################
###  Setting sample sizes  ###
##############################
## ASC fu
n_prob_asc <- 7063; n_sib_asc <- 2134

## SPARK Fu
n_prob_spk <- 5967; n_sib_spk <- 2615

## CNV
n_prob_cnv <- 5588+5733; n_sib_cnv <- 1667+2519

###################################
###  Loading gene count tables  ###
###################################
### SNV, denovo + inherited, ASC+GALA (sequenced at Broad Institute)
dg_asc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/merged-tables/oth_asc_fu_2023-11-06.txt')

### SNV, denovo + inherited, SPARK (sequenced at Regeneron)
dg_spk <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/merged-tables/oth_spark_fu_2023-11-06.txt')

dg_ddd <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/dataset-tables/dg_ddd_2024-02-08.txt')


## parameters ##
S <- sort(info$gene_id,index.return=TRUE)
  info <- info[S$ix,]
S <- sort(dg_asc$Gene_ID,index.return=TRUE)
  dg_asc <- dg_asc[S$ix,]
S <- sort(dg_spk$Gene_ID,index.return=TRUE)
  dg_spk <- dg_spk[S$ix,]
S <- sort(dg_ddd$Gene_ID,index.return=TRUE)
  dg_ddd <- dg_ddd[S$ix,]

##################################
###  calibrate mutation rates  ###
##################################
# BI-sequenced
cal_asc <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/calibration-tables/calibration-asc_2024-03-29.txt')

dg_asc$mut.ptv <- cal_asc$mut.ptv[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]
dg_asc$mut.misb <- cal_asc$mut.misb[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]
dg_asc$mut.misa <- cal_asc$mut.misa[match(dg_asc$Gene_ID, cal_asc$Gene_ID)]

# Regeneron-sequenced
cal_spk <- fread('/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/pop-max/calibration-tables/calibration-spk_2024-03-29.txt')

dg_spk$mut.ptv <- cal_spk$mut.ptv[match(dg_spk$Gene_ID, cal_spk$Gene_ID)]
dg_spk$mut.misb <- cal_spk$mut.misb[match(dg_spk$Gene_ID, cal_spk$Gene_ID)]
dg_spk$mut.misa <- cal_spk$mut.misa[match(dg_spk$Gene_ID, cal_spk$Gene_ID)]


#################################
### Loading CNV count tables  ###
#################################

### CNV, denovo, ASC+SSC+SPARK
cnv_dn <- fread('/sc/arion/projects/buxbaj01a/GALA/DATA/edited_data/TADA_August2023/non_AMR_CNV_counts_2023-08-25.txt')
cnv_dn <- cnv_dn[-which(cnv_dn$chr=="chr21" & (cnv_dn$end-cnv_dn$start)>30000000),] ### remove trisomy 21

###################################
###  Bayes Factor calculations  ###
###################################
beta.dn <- 0.2

###################
###  SNV/indel  ###
###################

###################
### de novo
BF_dn_ptv_asc <- BF_DN_SNV(count_case=dg_asc$dn.ptv, count_con=dg_asc$dn.ptv.sib, n_case=n_prob_asc, n_con=n_sib_asc, mut=dg_asc$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv_spk <- BF_DN_SNV(count_case=dg_spk$dn.ptv, count_con=dg_spk$dn.ptv.sib, n_case=n_prob_spk, n_con=n_sib_spk, mut=dg_spk$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

### DN misB
BF_dn_misB_asc <- BF_DN_SNV(count_case=dg_asc$dn.misb, count_con=dg_asc$dn.misb.sib, n_case=n_prob_asc, n_con=n_sib_asc, mut=dg_asc$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB_spk <- BF_DN_SNV(count_case=dg_spk$dn.misb, count_con=dg_spk$dn.misb.sib, n_case=n_prob_spk, n_con=n_sib_spk, mut=dg_spk$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

### DN misA
BF_dn_misA_asc <- BF_DN_SNV(count_case=dg_asc$dn.misa, count_con=dg_asc$dn.misa.sib, n_case=n_prob_asc, n_con=n_sib_asc, mut=dg_asc$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA_spk <- BF_DN_SNV(count_case=dg_spk$dn.misa, count_con=dg_spk$dn.misa.sib, n_case=n_prob_spk, n_con=n_sib_spk, mut=dg_spk$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk

###################
### inherited
BF_in_ptv <- BF_CC_SNV(count_case=dg_asc$t.ptv+dg_spk$t.ptv, count_con=dg_asc$u.ptv+dg_spk$u.ptv, n_case=n_prob_asc+n_prob_spk, n_con=n_prob_asc+n_prob_spk, mut=dg_asc$mut.ptv, prior=info$prior.in.ptv)
### IN misB
BF_in_misB <- BF_CC_SNV(count_case=dg_asc$t.misb+dg_spk$t.misb, count_con=dg_asc$u.misb+dg_spk$u.misb,  n_case=n_prob_asc+n_prob_spk, n_con=n_prob_asc+n_prob_spk, mut=dg_asc$mut.misb, prior=info$prior.in.misb)
### IN misA
BF_in_misA <- BF_CC_SNV(count_case=dg_asc$t.misa+dg_spk$t.misa, count_con=dg_asc$u.misa+dg_spk$u.misa,  n_case=n_prob_asc+n_prob_spk, n_con=n_prob_asc+n_prob_spk, mut=dg_asc$mut.misa, prior=info$prior.in.misa)

############
### CNV  ###
############
# size range, in number of constrained genes, to consider
cnv_size_range <- 1:8 

###################
### de novo
cnv <- processCNV(cnv_dn, loeuf_threshold=0.6,  info=info)
cnv_use <- cnv[which((is.na(cnv$gd_loci) | cnv$nahr==FALSE) & cnv$non_diploid_freq<.01 & cnv$num_genes %in% cnv_size_range & !is.na(cnv$Affected_Status)),]

  cnv_use$key <- paste(cnv_use$sample, cnv_use$chr, sep = '::')
  cnv_use <- cnv_use[cnv_use$key %notin% c('SP0036828::chr17')] # 240 --> 237
  cnv_use$key <- NULL

del_dup_adj <- table(cnv[,c("Affected_Status", "call")])
del_dup_adj <- (del_dup_adj[2,1]/del_dup_adj[1,1])/(del_dup_adj[2,2]/del_dup_adj[1,2])

del_use <- cnv_use[cnv_use$call=="DEL",]
mut.pred.del <- sapply(cnv_size_range, function(x){ max(1, sum(del_use$num_genes==x &del_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
del_use$mut <- mut.pred.del[del_use$num_genes]

BF_dn_del_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv)
BF_dn_del_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=del_use[del_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv)
BF_dn_del <- pmax(1, BF_dn_del_prob/BF_dn_del_sib)

dup_use <- cnv_use[cnv_use$call=="DUP",]
mut.pred.dup <- sapply(cnv_size_range, function(x){ max(1, sum(dup_use$num_genes==x &dup_use$Affected_Status==1))/n_sib_cnv/(length(info$prior.dn.ptv)-x)})
dup_use$mut <- mut.pred.dup[dup_use$num_genes]

BF_dn_dup_prob <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==2,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_prob_cnv, del_dup_adj=del_dup_adj)
BF_dn_dup_sib <- BF_DN_CNV(cnv_size_range=cnv_size_range, cnv_use=dup_use[dup_use$Affected_Status==1,], info=info, prior=info$prior.dn.ptv, beta.dn=beta.dn, n=n_sib_cnv, del_dup_adj=del_dup_adj)
BF_dn_dup <- pmax(1, BF_dn_dup_prob/BF_dn_dup_sib)


################################################
###  Accounting for genes with ACMG masking  ###
################################################
acmg_ind <- which(info$gene %in% acmg)

BF_dn_ptv_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.ptv, count_con=dg_spk$dn.ptv.sib, n_case= n_prob_spk, n_con=0, mut=dg_spk$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)
BF_dn_ptv_spk[acmg_ind] <- BF_dn_ptv_spk_acmg[acmg_ind]
BF_dn_ptv <- BF_dn_ptv_asc*BF_dn_ptv_spk

BF_dn_misB_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.misb, count_con=dg_spk$dn.misb.sib, n_case= n_prob_spk, n_con=0, mut=dg_spk$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)
BF_dn_misB_spk[acmg_ind] <- BF_dn_misB_spk_acmg[acmg_ind]
BF_dn_misB <- BF_dn_misB_asc*BF_dn_misB_spk

BF_dn_misA_spk_acmg <- BF_DN_SNV(count_case=dg_spk$dn.misa, count_con=dg_spk$dn.misa.sib, n_case= n_prob_spk, n_con=0, mut=dg_spk$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)
BF_dn_misA_spk[acmg_ind] <- BF_dn_misA_spk_acmg[acmg_ind]
BF_dn_misA <- BF_dn_misA_asc*BF_dn_misA_spk


########################################################################
### Integrating DDD - SNV/indel
########################################################################

### DN PTV
BF_dn_ptv_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.ptv, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.ptv, prior=info$prior.dn.ptv, beta.dn=beta.dn)

### DN misB
BF_dn_misB_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misb, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.misb, prior=info$prior.dn.misb, beta.dn=beta.dn)

### DN misA
BF_dn_misA_ddd <- BF_DN_SNV(count_case=dg_ddd$dn.misa, count_con=0, n_case=31058, n_con=0, mut=dg_ddd$mut.misa, prior=info$prior.dn.misa, beta.dn=beta.dn)


##################################
###  TADA outcome calculation  ###
##################################

### Aggregate Bayes Factors
# de novo only
BF_asd_dn <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk, 1),
                   pmax(BF_dn_misB_asc*BF_dn_misB_spk, 1),
                   pmax(BF_dn_misA_asc*BF_dn_misA_spk, 1),
                   pmax(BF_dn_del, 1), 
                   pmax(BF_dn_dup, 1))

# de novo + inherited
BF_asd_inh <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv, 1),
                    pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB, 1),
                    pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA, 1),
                    pmax(BF_dn_del, 1), 
                    pmax(BF_dn_dup, 1))

### Aggregate DD BF
BF_ddd <- cbind(pmax(1, BF_dn_ptv_ddd), pmax(BF_dn_misB_ddd, 1), pmax(BF_dn_misA_ddd, 1))

# Remove CNV evidence after aggregating 
### Aggregate NDD BF
BF_asd_ddd <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk*BF_in_ptv*BF_dn_ptv_ddd, 1),
                    pmax(BF_dn_misB_asc*BF_dn_misB_spk*BF_in_misB*BF_dn_misB_ddd, 1),
                    pmax(BF_dn_misA_asc*BF_dn_misA_spk*BF_in_misA*BF_dn_misA_ddd, 1),
                    pmax(BF_dn_del, 1), pmax(BF_dn_dup, 1))


# Only use del/dup if there is snv/indel evidence (Tot BF > 2) 
# ASD evidence
BF_asd_dn[apply(BF_asd_dn[,1:3], 1, max)<5,4:5] <- 1
BF_asd_inh[apply(BF_asd_inh[,1:3], 1, max)<5,4:5] <- 1
# ASD + DD evidence
BF_asd_ddd[apply(BF_asd_ddd[,1:3], 1, max)<5,4:5] <- 1

# Only use del/dup if there is snv/indel evidence (Tot BF > 2) 

### FDR calculation
# de novo only
qval_asd_dn <- Bayesian.FDR(apply(BF_asd_dn, 1, prod), pi0 = 1 - 0.05)
names(qval_asd_dn) <- info$gene_id

# de novo + inherited
qval_asd_inh <- Bayesian.FDR(apply(BF_asd_inh, 1, prod), pi0 = 1 - 0.05)
names(qval_asd_inh) <- info$gene_id

# DD
qval_ddd <- Bayesian.FDR(apply(BF_ddd, 1, prod), pi0 = 1 - 0.05)
names(qval_ddd) <- info$gene_id

# DD + OTHER
qval_asd_ddd <- Bayesian.FDR(apply(BF_asd_ddd, 1, prod), pi0=1-.05)
names(qval_asd_ddd) <- info$gene_id



### TADA-ASD significant genes at various thresholds
sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_dn<x))
# 56  87 130 169  de novo only for oth

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_inh<x))
# 59  92 140 187  de novo + inherited for oth

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_ddd<x))
# 309 378 477 559 DD with previous calibration

sapply(c(0.001, 0.01, 0.05, 0.1), function(x) sum(qval_asd_ddd<x))
# 358 457 608 737

#######################################

### back-transform FDR to p-values
# de novo only
ec_threshold <- .05/length(qval_asd_dn)
pval_asd_dn <- q2p(qval_asd_dn)
sum(pval_asd_dn <= ec_threshold) # 55

# de novo + inherited
ec_threshold <- .05/length(qval_asd_inh)
pval_asd_inh <- q2p(qval_asd_inh)
sum(qval_asd_inh <= ec_threshold) # 23

# DD
ec_threshold <- .05/length(qval_ddd)
pval_ddd <- q2p(qval_ddd)
sum(qval_ddd <= ec_threshold) # 212

# DD + GALA
ec_threshold <- .05/length(qval_asd_ddd)
pval_asd_ddd <- q2p(qval_asd_ddd)
sum(qval_asd_ddd <= ec_threshold) # 246

BF_table <- cbind(pmax(BF_dn_ptv_asc*BF_dn_ptv_spk, 1),
                  pmax(BF_dn_misB_asc*BF_dn_misB_spk, 1),
                  pmax(BF_dn_misA_asc*BF_dn_misA_spk, 1),
                  pmax(BF_dn_del, 1), 
                  pmax(BF_dn_dup, 1),
                  pmax(BF_in_ptv, 1), 
                  pmax(BF_in_misB, 1), 
                  pmax(BF_in_misA, 1))

df <- as.data.frame(BF_table)
df2 <- as.data.frame(BF_asd_inh)

names(df) <- c('BF.dn.ptv', 'BF.dn.misb', 'BF.dn.misa', 'BF.dn.del', 'BF.dn.dup', 
               'BF.in.ptv', 'BF.in.misb', 'BF.in.misa')

names(df2) <- c('BF.ptv', 'BF.misb', 'BF.misa', 'BF.del', 'BF.dup')


res <- data.table(gene = info$gene, gene_id = info$gene_id, pLI = info$pLI, LOEUF = info$LOEUF,
                  BF.dn.ptv = df$BF.dn.ptv, BF.dn.misb = df$BF.dn.misb, BF.dn.misa = df$BF.dn.misa, BF.dn.del = df$BF.dn.del, BF.dn.dup = df$BF.dn.dup,
                  BF.in.ptv = df$BF.in.ptv, BF.in.misB = df$BF.in.misb, BF.in.misA = df$BF.in.misa, 
                  qval_dn = qval_asd_dn, pval_dn = pval_asd_dn, qval_inh = qval_asd_inh, pval_inh = pval_asd_inh,
                  qval_gala_dd = qval_asd_ddd, pval_gala_dd = pval_asd_ddd, qval_dd = qval_ddd, pval_dd = pval_ddd,
                  l10PTV = log10(df2$BF.ptv), l10misb = log10(df2$BF.misb), l10misa = log10(df2$BF.misa), l10DEL = log10(df2$BF.del), l10DUP = log10(df2$BF.dup))

# write.table(res, '/sc/arion/projects/buxbaj01a/GALA/GALA_manuscript/DATA/edited_data/manuscript/2024/TADA_results/oth_bf_table_2024-04-17.txt', col.names = T, row.names = F, quote = F, sep = '\t')

