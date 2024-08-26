

library(data.table)
# setwd("~/Documents/Project/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/PRS/")
setwd("/Volumes/medpop_esp2/mesbah/projects/ch_progression/aric/gwas/PRS")

  # phenotype for incident CH
# aric_baseline_n_v05 <- fread("../../Display/aric_baseline_n_v05_N4187.pheno_ch_status.noHemeCA.correct_lipids.Jun3May2023.csv",
#                   header = T)
# aric_baseline_n_v05 <- fread("../../Display/aric_baseline_n_v05_N3730.pheno_ch_status.noHemeCA.correct_lipids.FinalDataset_4_glm.July132023.csv",
#                              header = T)
aric_baseline_n_v05 <- fread("../../pheno/aric_baseline_n_v05_N3730.pheno_ch_status.noHemeCA.correct_lipids.FinalDataset_4_glm.July132023.csv",
                             header = T)
# dbGaP linker
aric_dbgap_linker <- fread("../../pheno/phs000280.v7.pht001441.v7.p1.ARIC_Sample.MULTI.txt.gz", 
                           skip = 10, fill = TRUE, header = T)
aric_dbgap_linker <- aric_dbgap_linker[,c(1,2,4,5)]

aric_dbgap_linker <- subset(aric_dbgap_linker, !duplicated(aric_dbgap_linker$SAMPLE_ID))

table(aric_baseline_n_v05$GWAS_ID %in% aric_dbgap_linker$SAMPLE_ID, exclude = NULL)
# TRUE 
# 3730

# Phenotype for GWAS linker
# aric_baseline_n_v05_aa <- fread("../../aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv")
# aric_baseline_n_v05_ea <- fread("../../aric_baseline_n_v05_ea_PCA_GWAS.2023Feb20.tsv")
aric_baseline_n_v05_aa <- fread("../../aric_baseline_n_v05_aa_PCA_GWAS.2023Feb21.tsv")
aric_baseline_n_v05_ea <- fread("../../aric_baseline_n_v05_ea_PCA_GWAS.2023Feb21.tsv")

############
## PRS
############

############# CH
  ### CH EA
ch_ea <- fread("prs.ea_aric_rg22ch_beta.21SNV.sscore", 
                   header = T)
ch_ea <- merge(ch_ea[,c(2:5)], 
               aric_baseline_n_v05_ea[,c(2,3,47:56)], 
               by="IID")
ch_ea <- merge(ch_ea, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
  ### CH AA
ch_aa <- fread("prs.aa_aric_rg22ch_beta.21SNV.sscore", 
               header = T)
ch_aa <- merge(ch_aa[,c(2:5)], 
               aric_baseline_n_v05_aa[,c(2,3,47:56)], 
               by="IID")
ch_aa <- merge(ch_aa, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
################

####### Density plot of raw prs
library(ggplot2)
library(cowplot)
ch_aa <- ch_aa[, c(1:5)]
ch_aa$Ancestry <- "AA (n=637)"
ch_ea <- ch_ea[, c(1:5)]
ch_ea$Ancestry <- "EA (n=2376)"

ch_aa_ea <- as.data.frame(rbind(ch_aa, ch_ea))

str(ch_aa_ea)

# Use semi-transparent fill
pdf("Fig4a.DistributionofprsCH.pdf", 
    width = 7, height= 7)

ggplot(ch_aa_ea, aes(x=SCORE1_AVG, fill=Ancestry)) +
  geom_density(alpha=0.4) + xlab("Prevalent CH PRS")

dev.off()
######

################ DNMT3A
### DNMT3A EA
dnmt3a_ea <- fread("prs.ea_aric_rg22dnmt3a_beta.22SNV.sscore", 
               header = T)
dnmt3a_ea <- merge(dnmt3a_ea[,c(2:5)], 
               aric_baseline_n_v05_ea[,c(2,3,47:56)], 
               by="IID")
dnmt3a_ea <- merge(dnmt3a_ea, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
### DNMT3A AA
dnmt3a_aa <- fread("prs.aa_aric_rg22dnmt3a_beta.22SNV.sscore", 
               header = T)
dnmt3a_aa <- merge(dnmt3a_aa[,c(2:5)], 
               aric_baseline_n_v05_aa[,c(2,3,47:56)], 
               by="IID")
dnmt3a_aa <- merge(dnmt3a_aa, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
###############

############### TET2
### TET2 EA
tet2_ea <- fread("prs.ea_aric_rg22tet2_beta.6SNV.sscore", 
               header = T)
tet2_ea <- merge(tet2_ea[,c(2:5)], 
               aric_baseline_n_v05_ea[,c(2,3,47:56)], 
               by="IID")
tet2_ea <- merge(tet2_ea, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
### TET2 AA
tet2_aa <- fread("prs.aa_aric_rg22tet2_beta.6SNV.sscore", 
               header = T)
tet2_aa <- merge(tet2_aa[,c(2:5)], 
               aric_baseline_n_v05_aa[,c(2,3,47:56)], 
               by="IID")
tet2_aa <- merge(tet2_aa, 
               aric_baseline_n_v05, 
               by="GWAS_ID")
#####################

##### Save 
# save.image("/Volumes/medpop_esp2/mesbah/projects/ch_progression/aric/gwas/PRS/PRS.Prev_CH.n_pheno.2023Jul14.rda")
##### PRS and incident CH
load("/Volumes/medpop_esp2/mesbah/projects/ch_progression/aric/gwas/PRS/PRS.Prev_CH.n_pheno.2023Jul14.rda")


#####
##############################
#### Forest plot
##############################
library(data.table) # version 1.14.6
library(meta) # version 6.2-1
library(grid) # version 4.2.2
library(scales) # version 1.2.1
##############################
################### Inverse-variance weighted meta-analysis ####
# invMeta2 <- function(b1, b2, bSE1, bSE2){
#   w1 <- 1/(bSE1^2)
#   w2 <- 1/(bSE2^2)
#   
#   se <- sqrt(1/sum(w1,w2))
#   beta <- sum(b1*w1, b2*w2)/sum(w1, w2)
#   
#   Z <- beta/se
#   
#   pval <- 2 * pnorm(abs(Z), lower.tail = FALSE)
#   
#   # pnorm(abs(d3$metaBeta)/d3$metaSE, lower.tail = FALSE)
#   return(list(beta, se, Z, pval))
# }

source("/Volumes/medpop_esp2/mesbah/tools/longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/ivw_metaAnalysis.R")


# set working directory
setwd("/Volumes/medpop_esp2/mesbah/projects/ch_progression/aric/gwas/PRS")
# PRS summary
prs_dat <- fread("ch_prs_std.final_glm.multivariable.incident_ch.2023Jul15.csv", 
                 header = T)

table(prs_dat$Dataset)
table(prs_dat$Exposure)
table(prs_dat$Outcome)
prs_dat.vaf2 <- subset(prs_dat, 
                       prs_dat$Dataset %in% 
                         c("AA", "EA"))

prs_dat.vaf10 <- subset(prs_dat, 
                       prs_dat$Dataset %in% 
                         c("AA_vaf10", "EA_vaf10")) 

### ivw meta-analysis ###
meta_prs.vaf2 <- as.data.frame(matrix(NA, 
                                 nrow = 4, 
                                 ncol = ncol(prs_dat.vaf2)))
names(meta_prs.vaf2) <- names(prs_dat.vaf2)
meta_prs.vaf2$Outcome <- unique(prs_dat.vaf2$Outcome)
meta_prs.vaf2$Dataset <- "Overall (VAF>=2%)"
meta_prs.vaf2$Exposure <- "CH_PRS"

# invMeta2(b1 = adj.glm_dat_v1$Beta[adj.glm_dat_v1$Outcome=="incident_CH"][1],
#          b2 = adj.glm_dat_v1$Beta[adj.glm_dat_v1$Outcome=="incident_CH"][2],
#          bSE1 = adj.glm_dat_v1$SE[adj.glm_dat_v1$Outcome=="incident_CH"][1], 
#          bSE2 = adj.glm_dat_v1$SE[adj.glm_dat_v1$Outcome=="incident_CH"][2])


for(i in 1:nrow(meta_prs.vaf2)){
  tmp <- ivw_meta(b = prs_dat.vaf2$Beta[prs_dat.vaf2$Outcome==meta_prs.vaf2$Outcome[i]], 
                  bSE = prs_dat.vaf2$SE[prs_dat.vaf2$Outcome==meta_prs.vaf2$Outcome[i]])
  meta_prs.vaf2$Beta[i] <- tmp[[1]]
  meta_prs.vaf2$SE[i] <- tmp[[2]]
  meta_prs.vaf2$`t-stat`[i] <- tmp[[3]]
  meta_prs.vaf2$P[i] <- tmp[[4]]
}

prs_dat_n_meta.vaf2 <- as.data.frame(rbind(prs_dat.vaf2,
                                            meta_prs.vaf2))

### VAF>=10%
meta_prs.vaf10 <- as.data.frame(matrix(NA, 
                                      nrow = 4, 
                                      ncol = ncol(prs_dat.vaf10)))
names(meta_prs.vaf10) <- names(prs_dat.vaf10)
meta_prs.vaf10$Outcome <- unique(prs_dat.vaf10$Outcome)
meta_prs.vaf10$Dataset <- "Overall (VAF>=10%)"
meta_prs.vaf10$Exposure <- "CH_PRS"

for(i in 1:nrow(meta_prs.vaf10)){
  tmp <- ivw_meta(b = prs_dat.vaf10$Beta[prs_dat.vaf10$Outcome==meta_prs.vaf10$Outcome[i]], 
                  bSE = prs_dat.vaf10$SE[prs_dat.vaf10$Outcome==meta_prs.vaf10$Outcome[i]])
  meta_prs.vaf10$Beta[i] <- tmp[[1]]
  meta_prs.vaf10$SE[i] <- tmp[[2]]
  meta_prs.vaf10$`t-stat`[i] <- tmp[[3]]
  meta_prs.vaf10$P[i] <- tmp[[4]]
}


prs_dat_n_meta.vaf10 <- as.data.frame(rbind(prs_dat.vaf10,
                                     meta_prs.vaf10))

###
## Outcome
prs_dat_n_meta.vaf2$OUTCOME[prs_dat_n_meta.vaf2$Outcome=="incident_CH"] <- "Overall CH"
prs_dat_n_meta.vaf2$OUTCOME[prs_dat_n_meta.vaf2$Outcome=="incident_DNMT3A"] <- "DNMT3A"
prs_dat_n_meta.vaf2$OUTCOME[prs_dat_n_meta.vaf2$Outcome=="incident_TET2"] <- "TET2"
prs_dat_n_meta.vaf2$OUTCOME[prs_dat_n_meta.vaf2$Outcome=="incident_ASXL1"] <- "ASXL1"

prs_dat_n_meta.vaf2$Exp[prs_dat_n_meta.vaf2$Dataset=="AA"] <- "ARIC AA"
prs_dat_n_meta.vaf2$Exp[prs_dat_n_meta.vaf2$Dataset=="EA"] <- "ARIC EA"
prs_dat_n_meta.vaf2$Exp[prs_dat_n_meta.vaf2$Exp=="CH_PRS"] <- "Overall"
table(prs_dat_n_meta.vaf2$OUTCOME)
table(prs_dat_n_meta.vaf2$Exp)
####
prs_dat_n_meta.vaf2$color <- "gray"
prs_dat_n_meta.vaf2$color[prs_dat_n_meta.vaf2$Exp=="Overall"] <- "red"
# box shape
prs_dat_n_meta.vaf2$Types <- "square"
prs_dat_n_meta.vaf2$Types[prs_dat_n_meta.vaf2$Exp=="Overall"] <- "diamond"


##

##################

# ## Exposures
# prs_dat$Exposure[prs_dat$Exposure=="AA_CH_PRS"] <- "Prevalent CH PRS"
# 
# prs_dat$Exposure[prs_dat$Exposure=="EA_CH_PRS"] <- "Prevalent CH PRS"
# 
# table(prs_dat$Exposure)



# ## 20 independent test at 5%; P< 0.05/20 = 0.0025
# # cat("P threshold< 0.0025")
# # 0.05/20 = 0.0025 = "2.5E-03"
# # 0.05/15 = 0.0033 = "3.3E-03"
# cat("P threshold< ",round(0.05/15,4))
# prs_dat$sig <- ifelse(prs_dat$P<0.0033, "***","")
# table(prs_dat$sig)

# format 
prs_dat_n_meta.vaf2$P_val <- formatC(x = prs_dat_n_meta.vaf2$P, 
                         digits = 1,
                         format = "E")

# OR
prs_dat_n_meta.vaf2$OR <- formatC(round(exp(prs_dat_n_meta.vaf2$Beta),2), 
                      digits = 2, format = "f")

prs_dat_n_meta.vaf2$lSE <- ( prs_dat_n_meta.vaf2$Beta - 1.96 * prs_dat_n_meta.vaf2$SE)
prs_dat_n_meta.vaf2$uSE <- ( prs_dat_n_meta.vaf2$Beta + 1.96 * prs_dat_n_meta.vaf2$SE)

# 95% CI
prs_dat_n_meta.vaf2$CI95 <- paste0("[",formatC(round(exp( prs_dat_n_meta.vaf2$Beta - 1.96 * prs_dat_n_meta.vaf2$SE),2), 
                                               digits = 2, format = "f"),
                       ", ",
                       formatC(round(exp( prs_dat_n_meta.vaf2$Beta + 1.96 * prs_dat_n_meta.vaf2$SE),2), 
                               digits = 2, format = "f"), 
                       "]")
(prs_dat_n_meta.vaf2)

# 
# ## Adjusted
# adj.glm_dat_v1 <- subset(prs_dat, 
#                          prs_dat$Dataset %in% 
#                            c("AA", "EA"))
## Sort outcome
prs_dat_n_meta.vaf2$OUTCOME <- ordered(prs_dat_n_meta.vaf2$OUTCOME, 
                                  levels = c("Overall CH", 
                                             "DNMT3A", "TET2", 
                                             "ASXL1"))
## exposure
prs_dat_n_meta.vaf2$Exp <- ordered(prs_dat_n_meta.vaf2$Exp,
                              levels = c("ARIC AA", 
                                         "ARIC EA", 
                                         "Overall"))

## cox_dat_v1: CH, DNMT3A, TET2
## All adjusted
# df_plot_adj_v1 <- subset(glm_dat_v1, glm_dat_v1$Dataset=="Adjusted")
## Format metagen
b_adj_v1 <- metagen(TE = Beta,
                    lower = lSE,
                    upper = uSE,
                    studlab = Exp,
                    subgroup=OUTCOME,
                    data=prs_dat_n_meta.vaf2,
                    sm="OR", 
                    common = F,
                    overall = F, 
                    overall.hetstat = F)

### adjusted CH, DNMT3A, TET2, ASXL1
pdf("prsCH.final_glm.Forest_incidentCH.2023Jul17.pdf", 
    width = 7, height= 7)
forest(x = b_adj_v1,
       sortvar = prs_dat_n_meta.vaf2$Exp,
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="common",  
       level=0.95, 
       xlim=c(0.5, 3), 
       smlab="Effect of prevalent CH PRS\non Incident CH\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio (95% CI) for incident CH\nper SD increase in prevalent CH PRS", 
       squaresize=0.6, 
       col.subgroup="black",
       col.square = prs_dat_n_meta.vaf2$color,
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab"), 
       leftlabs = c("                     "),
       rightcols=c("OR","CI95","P_val"),
       rightlabs=c("OR","95% CI","P"),
       #rightcols=NULL, 
       #rightlabs=NULL,
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)
dev.off()


######### library
# library(pdftools)
# library(png)
# library(patchwork)
library(magick)
prs_png <- image_read("prsCH.final_glm.Forest_incidentCH.2023Jul17.png")
gwas_png <- image_read("../GWAS.Forest_incidentCH.v2.28May23.italic.png")


# Resize the images to the same dimensions (if needed)
prs_png <- image_scale(prs_png, "1400x1400")
gwas_png <- image_scale(gwas_png, "2000x1600")

# Combine the two images side by side
combined_image <- image_append(c(prs_png , gwas_png), 
                               stack = F)

# Save the combined image as a new PNG file
image_write(combined_image, 
            "fig4ab.combined_image.side.2023Jul17.png")

