
###### GWAS Overlap
setwd("/ch_progression/aric/gwas/gwas_current/meta_aric/")
library(data.table)
library(dplyr)
rg22_ch_dnmt_tet <- fread("../../rg2022/rg22_ch_dnmt_tet.hg19.csv", header=T)
counts.var <- as.data.frame(table(rg22_ch_dnmt_tet$varID_hg19))
var1.rg22_ch_dnmt_tet <- rg22_ch_dnmt_tet[rg22_ch_dnmt_tet$varID_hg19 %in% counts.var$Var1[counts.var$Freq==1], ]

# var2.rg22_ch_dnmt_tet <- rg22_ch_dnmt_tet[rg22_ch_dnmt_tet$varID_hg19 %in% counts.var$Var1[counts.var$Freq>1], 
var2.rg22_ch_dnmt_tet <- rg22_ch_dnmt_tet[, c(1,2,3,5,
                                            6,9,10,
                                            12,15)] %>% 
  group_by(varID_hg19) %>% 
  summarise(varID_hg38=unique(varID),
            chrom=unique(chr_hg19),
            POS_hg19=unique(POS_hg19),
            Traits = paste(Trait, collapse = "|"),
            Effect_UKB=max(Effect_UKB),
            Pval_UKB=min(Pval_UKB),
            AAF_UKB=max(AAF_UKB),
            nearestGene=unique(nearestGene))


  # CH
ch <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_CH.wao.tsv.gz |  cut -f1-38  ", header = F)
names(ch) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

ch$varID_aric <- paste(ch$chr_1, ch$POS, 
                          ch$REF, ch$ALT, sep=":")
table(ch$varID_aric==ch$varID_hg19)
# 8

ch_growth <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_CH_or_growingClones.wao.tsv.gz |  cut -f1-38  ", header = F)
names(ch_growth) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

ch_growth$varID_aric <- paste(ch_growth$chr_1, ch_growth$POS, 
                              ch_growth$REF, ch_growth$ALT, sep=":")
table(ch_growth$varID_aric==ch_growth$varID_hg19)
# 13

dta <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_DTA.wao.tsv.gz |  cut -f1-38  ", header = F)
names(dta) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

dta$varID_aric <- paste(dta$chr_1, dta$POS, 
                          dta$REF, dta$ALT, sep=":")
table(dta$varID_aric==dta$varID_hg19)
# 24

sf <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_SF.wao.tsv.gz |  cut -f1-38  ", header = F)
names(sf) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

sf$varID_aric <- paste(sf$chr_1, sf$POS, 
                       sf$REF, sf$ALT, sep=":")
table(sf$varID_aric==sf$varID_hg19)
# 6

ddr <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_DDR.wao.tsv.gz |  cut -f1-38  ", header = F)
names(ddr) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

ddr$varID_aric <- paste(ddr$chr_1, ddr$POS, 
                        ddr$REF, ddr$ALT, sep=":")
table(ddr$varID_aric==ddr$varID_hg19)
# 2


dnmt <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_DNMT3A.wao.tsv.gz |  cut -f1-38  ", header = F)

names(dnmt) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                "chr_1","POS_1","POS_2","chr_2", "POS", 
                "REF", "ALT",	"reference_allele",	"other_allele",
                "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                "z",	"P",	"log10_P", "q_statistic",	
                "q_p",	"i2",	"n_studies",	"n_samples",
                "effects")

dnmt$varID_aric <- paste(dnmt$chr_1, dnmt$POS, 
                         dnmt$REF, dnmt$ALT, sep=":")
table(dnmt$varID_aric==dnmt$varID_hg19)
# 23

tet <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_TET2.wao.tsv.gz |  cut -f1-38  ", header = F)
names(tet) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                "chr_1","POS_1","POS_2","chr_2", "POS", 
                "REF", "ALT",	"reference_allele",	"other_allele",
                "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                "z",	"P",	"log10_P", "q_statistic",	
                "q_p",	"i2",	"n_studies",	"n_samples",
                "effects")

tet$varID_aric <- paste(tet$chr_1, tet$POS, 
                        tet$REF, tet$ALT, sep=":")
table(tet$varID_aric==tet$varID_hg19)
# 15

  # ASLX1
asxl1 <- fread(cmd="gzcat overlap.p10.meta_aric_N3015.chr1_22.ea2378_aa637.incident_ASXL1.wao.tsv.gz |  cut -f1-38  ", header = F)
names(asxl1) <- c("chr_hg19",	"Start_500kb",	"End_500kb",	
                  "varID",	"POS_hg19",	"locus",	"Trait","Effect_UKB",
                  "LCI_Effect_UKB",	"UCI_Effect_UKB",	"Pval_UKB",	
                  "AAF_UKB",	"CHROM",	"varID_hg19",	"nearestGene",
                  "chr_1","POS_1","POS_2","chr_2", "POS", 
                  "REF", "ALT",	"reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")

asxl1$varID_aric <- paste(asxl1$chr_1, asxl1$POS, 
                          asxl1$REF, asxl1$ALT, sep=":")
table(asxl1$varID_aric==asxl1$varID_hg19)
asxl1[(asxl1$varID_aric==asxl1$varID_hg19),]

##### named variants
## while read pheno; do zgrep -wf <(awk 'NR>1{print $14}' /ch_progression/aric/gwas/rg2022/rg22_ch_dnmt_tet.hg19.tsv | sort -V | uniq) /ch_progression/aric/gwas/gwas_current/meta_aric/meta_aric_N3015.chr1_22.ea2378_aa637.incident_${pheno}.out.gz > /ch_progression/aric/gwas/gwas_current/meta_aric/rg22_variants.${pheno}.tsv; done < <(echo -e "CH\nCH_or_growingClones\nDTA\nSF\nDDR\nDNMT3A\nTET2\nASXL1" ) &
rg_ch <- fread("rg22_variants.CH.tsv", header = F)
names(rg_ch) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_ch <- merge(rg_ch, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_ch_pls <- fread("rg22_variants.CH_or_growingClones.tsv", header = F)
names(rg_ch_pls) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_ch_pls <- merge(rg_ch_pls, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_dta <- fread("rg22_variants.DTA.tsv", header = F)
names(rg_dta) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_dta <- merge(rg_dta, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_sf <- fread("rg22_variants.SF.tsv", header = F)
names(rg_sf) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_sf <- merge(rg_sf, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_ddr <- fread("rg22_variants.DDR.tsv", header = F)
names(rg_ddr) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_ddr <- merge(rg_ddr, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_dnmt <- fread("rg22_variants.DNMT3A.tsv", header = F)
names(rg_dnmt) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_dnmt <- merge(rg_dnmt, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_tet <- fread("rg22_variants.TET2.tsv", header = F)
names(rg_tet) <- c("varID","reference_allele",	"other_allele",
                  "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                  "z",	"P",	"log10_P", "q_statistic",	
                  "q_p",	"i2",	"n_studies",	"n_samples",
                  "effects")
rg_tet <- merge(rg_tet, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")

rg_asxl <- fread("rg22_variants.ASXL1.tsv", header = F)
names(rg_asxl) <- c("varID","reference_allele",	"other_allele",
                   "eaf",	"beta",	"se",	"beta_95L",	"beta_95U",	
                   "z",	"P",	"log10_P", "q_statistic",	
                   "q_p",	"i2",	"n_studies",	"n_samples",
                   "effects")

rg_asxl <- merge(rg_asxl, rg22_ch_dnmt_tet, 
               by.x="varID", by.y="varID_hg19")
##### Write excel file
require(openxlsx)

# Lead variants from RG22 GWAS
list_of_datasets <- list("CH"=rg_ch,"CH_pls_grt" =rg_ch_pls, 
                         "DTA"=rg_dta, "SF"=rg_sf, "DDR"=rg_ddr, 
                         "DNMT3A"=rg_dnmt, "TET2"=rg_tet,"ASXL1"=rg_asxl)
write.xlsx(list_of_datasets, file = "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableS.rg22_lead_vars.xlsx")


# all p<=0.1 within +-500k of RG22 GWAS loci
all_p10.ist_of_datasets <- list("CH"=ch,"CH_pls_grt" =ch_growth, 
                         "DTA"=dta, "SF"=sf, "DDR"=ddr, 
                         "DNMT3A"=dnmt, "TET2"=tet,"ASXL1"=asxl1)
write.xlsx(all_p10.ist_of_datasets, file = "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableS.aric__p10.rg22_lead_pls500kb.xlsx")

## get LD 
# https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html
# need RS id 
# install.packages("LDlinkR")
library(LDlinkR)
LDproxy(snp = "chr14:96180685",
        pop = "CEU", 
        r2d = "r2",
        token = "94afd01dbc34",
        file = "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TET2.TCL1A.LD.tsv", 
        genome_build = "grch37")

tcl1a_tet <- LDproxy(snp = "chr14:96180685",
        pop = "CEU", 
        r2d = "r2",
        token = "94afd01dbc34",
        genome_build = "grch37")
## pfile: /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --pfile /broad/hptmp/mesbah/dataset/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr14.1kg_imp --ld rs3087688 rs11846938

## /medpop/esp2/mesbah/tools/plink2_linux_x86_64_20200831 --pgen /broad/hptmp/mesbah/dataset/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr14.1kg_imp.pgen --psam /broad/hptmp/mesbah/dataset/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr14.1kg_imp.psam --pvar /broad/hptmp/mesbah/dataset/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr14.1kg_imp.pvar  --ld rs3087688 rs11846938
## 
###########
setwd("/ch_progression/aric/gwas/gwas_current/meta_aric/")
rg22.chip_hg38tohg19 <- fread("../../rg2022/ch_rg22.hg19.csv", 
                              header = T, sep = ",")
names(rg22.chip_hg38tohg19)
rg22.chip_hg38tohg19 <- rg22.chip_hg38tohg19[, c(1:4,9,11,12,13,14,15,48:51,38)]

rg22.dnmt3a_hg38tohg19 <- fread("../../rg2022/dnmt3a_rg22.hg19.csv", 
                                header = T, sep = ",")
names(rg22.dnmt3a_hg38tohg19) 
rg22.dnmt3a_hg38tohg19 <- rg22.dnmt3a_hg38tohg19[, c(1:4,9,11,12,13,14,15,45:48,38)]
names(rg22.dnmt3a_hg38tohg19) <-  names(rg22.chip_hg38tohg19)

rg22.tet2_hg38tohg19 <- fread("../../rg2022/tet2_rg22.hg19.csv", 
                              header = T, sep = ",")
rg22.tet2_hg38tohg19$locus <- c("11-1","14-1","14-1","17-1","3-1","3-2","5-1","5-1")
names(rg22.tet2_hg38tohg19)
rg22.tet2_hg38tohg19 <- rg22.tet2_hg38tohg19[, c(1:3,48,8,10:14,44:47,37)]
names(rg22.tet2_hg38tohg19) <-  names(rg22.chip_hg38tohg19)


rg22_ch_dnmt_tet <- as.data.frame(rbind(rg22.chip_hg38tohg19, rg22.dnmt3a_hg38tohg19, rg22.tet2_hg38tohg19))

# fwrite(rg22_ch_dnmt_tet, "../../rg2022/rg22_ch_dnmt_tet.hg19.csv", 
    #  col.names = T, row.names = F, sep = ",")

############
# CH phenotypes for GWAS
setwd("/ch_progression/aric/gwas/rg2022/")
# load prev. CH GWAS loci 
chip_hg38tohg19 <- fread("ch_rg22.hg19.csv", header=T)
str(chip_hg38tohg19)  

dnmt3a_hg38tohg19 <- fread("dnmt3a_rg22.hg19.csv", header=T)
str(dnmt3a_hg38tohg19)

tet2_hg38tohg19 <- fread("tet2_rg22.hg19.csv", header=T)
str(tet2_hg38tohg19)

## Load phenotypes
aric_baseline_n_v05 <- fread("../../pheno/aric_baseline_n_v05_N4189.pheno_ch_status_trajectory.23Mar2023.csv", header=T)
aric_baseline_n_v05$dAge <- aric_baseline_n_v05$Age - aric_baseline_n_v05$age_base

# dbGaP linker
aric_dbgap_linker <- fread("../../pheno/phs000280.v7.pht001441.v7.p1.ARIC_Sample.MULTI.txt.gz", 
                           skip = 10, fill = TRUE, header = T)

aric_dbgap_linker <- aric_dbgap_linker[,c(1,2,4,5)]

aric_dbgap_linker <- subset(aric_dbgap_linker, !duplicated(aric_dbgap_linker$SAMPLE_ID))

table(aric_baseline_n_v05$GWAS_ID %in% aric_dbgap_linker$SAMPLE_ID, exclude = NULL)


# link dbGaP ids with GWASID
aric_baseline_n_v05 <- merge(aric_dbgap_linker, 
                             aric_baseline_n_v05,
                             by.x = "SAMPLE_ID",
                             by.y = "GWAS_ID")

aric_baseline_n_v05$SUBJECT_ID <- as.numeric(aric_baseline_n_v05$SUBJECT_ID)

head(aric_baseline_n_v05)

## PCA data
pca_ea <- fread("../../pca/ARIC_EurAm_chr1_22.aray_snps.eigenvec", header = T)
# pca_ea$IID <- as.character(pca_ea$IID)
pca_aa <- fread("../../pca/ARIC_AfrAm_chr1_22.aray_snps.eigenvec", header = T)
# pca_ea$IID <- as.character(pca_ea$IID)
table(aric_baseline_n_v05$SUBJECT_ID %in% pca_ea$IID)

# EA: 2692
aric_baseline_n_v05_ea <- merge(aric_baseline_n_v05, pca_ea, 
                                by.x="SUBJECT_ID", by.y="IID")

nrow(aric_baseline_n_v05_ea)

# EA: 2692
aric_baseline_n_v05_aa <- merge(aric_baseline_n_v05, pca_aa, 
                                by.x="SUBJECT_ID", by.y="IID")

nrow(aric_baseline_n_v05_aa)
  # EA
names(aric_baseline_n_v05_ea[, c(1,1,2,6,3:5,55,9:21,109, 61,107,108,112:121,85,86,87,89)])
aric_baseline_n_v05_ea_gwas <- aric_baseline_n_v05_ea[, c(1,1,2,6,3:5,55,9:21,109, 61,107,108,112:121,85,86,87,89)]
names(aric_baseline_n_v05_ea_gwas) <- c("FID","IID","GWAS_ID", names(aric_baseline_n_v05_ea_gwas)[4:39])
str(aric_baseline_n_v05_ea_gwas)

  # AA
names(aric_baseline_n_v05_aa)
names(aric_baseline_n_v05_aa[, c(1,1,2,6,3:5,55,9:21,109, 61,107,108,112:121,85,86,87,89)])
aric_baseline_n_v05_aa_gwas <- aric_baseline_n_v05_aa[, c(1,1,2,6,3:5,55,9:21,109, 61,107,108,112:121,85,86,87,89)]
names(aric_baseline_n_v05_aa_gwas) <- c("FID","IID","GWAS_ID", names(aric_baseline_n_v05_aa_gwas)[4:39])
str(aric_baseline_n_v05_aa_gwas)

## GWAS_ID==SAMPLE_ID
fwrite(aric_baseline_n_v05_ea, "../aric_baseline_n_v05_ea_PCA.2023Apr3.csv", 
       row.names = F, col.names = T, sep=",", na="NA")
fwrite(aric_baseline_n_v05_ea_gwas, "../aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.tsv",
       row.names = F, col.names = T, sep="\t", na="NA", quote = F)
# 
fwrite(aric_baseline_n_v05_aa, "../aric_baseline_n_v05_aa_PCA.2023Apr3.csv", 
       row.names = F, col.names = T, sep=",", na="NA")
fwrite(aric_baseline_n_v05_aa_gwas, "../aric_baseline_n_v05_aa_PCA_GWAS.2023Apr3.tsv",
       row.names = F, col.names = T, sep="\t", na="NA", quote = F)
##############################################
##### RG 2022 GWAS 
library(readxl)
library(data.table)
ch_rg22 <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/41586_2022_5448_MOESM5_ESM.xlsx", 
                      sheet = 3)
dnmt3a_rg22 <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/41586_2022_5448_MOESM5_ESM.xlsx", 
                      sheet = 12)
tet2_rg22 <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/41586_2022_5448_MOESM5_ESM.xlsx", 
                      sheet = 15)
# BiocManager::install("liftOver")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# BiocManager::install("rtracklayer")
library(liftOver)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# tx_hg19 <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# unzipped chain file needed
# Load Chain file
hg38_to_19_chain = import.chain( "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/hg38ToHg19.over.chain")

ch_rg22$CHROM <- paste0("chr",ch_rg22$Chr)
hg38.chip_df <- ch_rg22[,c(46, 4, 4,2)]
names(hg38.chip_df) <- c("chrom", "start", "end","varID")
gr.chip.hg38 <- as(hg38.chip_df, "GenomicRanges")

dnmt3a_rg22$CHROM <- paste0("chr",dnmt3a_rg22$chr)
hg38.dnmt3a_df <- dnmt3a_rg22[,c(43, 4, 4,2)]
names(hg38.dnmt3a_df) <- c("chrom", "start", "end","varID")
gr.dnmt3a.hg38 <- as(hg38.dnmt3a_df, "GenomicRanges")

tet2_rg22$CHROM <- paste0("chr",tet2_rg22$chr)
hg38.tet2_df <- tet2_rg22[,c(42, 3, 3,1)]
names(hg38.tet2_df) <- c("chrom", "start", "end","varID")
gr.tet2.hg38 <- as(hg38.tet2_df, "GenomicRanges")

chip_hg38tohg19 <- liftOver(gr.chip.hg38, hg38_to_19_chain)
chip_hg38tohg19 <- as.data.frame(chip_hg38tohg19)
chip_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",replacement = "", 
                           x = as.character(chip_hg38tohg19$seqnames))
chip_hg38tohg19$POS_hg19 <- chip_hg38tohg19$start
chip_hg38tohg19 <- merge(chip_hg38tohg19[,c(8,9,10)], ch_rg22, 
                         by.x="varID",by.y="Name")
chip_hg38tohg19$varID_hg19 <- paste(chip_hg38tohg19$chr_hg19, 
                                    chip_hg38tohg19$POS_hg19, 
                                    chip_hg38tohg19$Ref, 
                                    chip_hg38tohg19$Alt, sep=":")
## DNMT3A
dnmt3a_rg22$CHROM <- paste0("chr",dnmt3a_rg22$chr)
hg38.dnmt3a_df <- dnmt3a_rg22[,c(43, 4, 4,2)]
names(hg38.dnmt3a_df) <- c("chrom", "start", "end","varID")
gr.dnmt3a.hg38 <- as(hg38.dnmt3a_df, "GenomicRanges")

dnmt3a_hg38tohg19 <- liftOver(gr.dnmt3a.hg38, hg38_to_19_chain)
dnmt3a_hg38tohg19 <- as.data.frame(dnmt3a_hg38tohg19)
dnmt3a_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",replacement = "", 
                                 x = as.character(dnmt3a_hg38tohg19$seqnames))
dnmt3a_hg38tohg19$POS_hg19 <- dnmt3a_hg38tohg19$start
dnmt3a_hg38tohg19 <- merge(dnmt3a_hg38tohg19[,c(8,9,10)], dnmt3a_rg22, 
                         by.x="varID",by.y="name")
dnmt3a_hg38tohg19$varID_hg19 <- paste(dnmt3a_hg38tohg19$chr_hg19, 
                                      dnmt3a_hg38tohg19$POS_hg19, 
                                      dnmt3a_hg38tohg19$ref, 
                                      dnmt3a_hg38tohg19$alt, 
                                      sep=":")

## TET2
tet2_rg22$CHROM <- paste0("chr",tet2_rg22$chr)
hg38.tet2_df <- tet2_rg22[,c(42, 3, 3,1)]
names(hg38.tet2_df) <- c("chrom", "start", "end","varID")
gr.tet2.hg38 <- as(hg38.tet2_df, "GenomicRanges")

tet2_hg38tohg19 <- liftOver(gr.tet2.hg38, hg38_to_19_chain)
tet2_hg38tohg19 <- as.data.frame(tet2_hg38tohg19)
tet2_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",replacement = "", 
                                 x = as.character(tet2_hg38tohg19$seqnames))
tet2_hg38tohg19$POS_hg19 <- tet2_hg38tohg19$start

tet2_hg38tohg19 <- merge(tet2_hg38tohg19[,c(8,9,10)], tet2_rg22, 
                         by.x="varID",by.y="name")
tet2_hg38tohg19$varID_hg19 <- paste(tet2_hg38tohg19$chr_hg19, 
                                    tet2_hg38tohg19$POS_hg19, 
                                    tet2_hg38tohg19$ref, 
                                    tet2_hg38tohg19$alt, sep=":")

chip_hg38tohg19$Start_500kb <- chip_hg38tohg19$POS_hg19-500000
chip_hg38tohg19$End_500kb <- chip_hg38tohg19$POS_hg19+500000
dnmt3a_hg38tohg19$Start_500kb <- dnmt3a_hg38tohg19$POS_hg19-500000
dnmt3a_hg38tohg19$End_500kb <- dnmt3a_hg38tohg19$POS_hg19+500000
tet2_hg38tohg19$Start_500kb <- tet2_hg38tohg19$POS_hg19-500000
tet2_hg38tohg19$End_500kb <- tet2_hg38tohg19$POS_hg19+500000

fwrite(chip_hg38tohg19, "/Meta_GWAS/2023_gwas/rg2022/ch_rg22.hg19.csv", 
       row.names = F, col.names = T, sep = ",")
fwrite(dnmt3a_hg38tohg19, "/Meta_GWAS/2023_gwas/rg2022/dnmt3a_rg22.hg19.csv", 
       row.names = F, col.names = T, sep = ",")
fwrite(tet2_hg38tohg19, "/Meta_GWAS/2023_gwas/rg2022/tet2_rg22.hg19.csv", 
       row.names = F, col.names = T, sep = ",")

save.image("/Meta_GWAS/2023_gwas/rg2022/rg22_gwas_loci.hg19.rda")
######################### 

### GWAS catalog summary 
ch_regen <- fread("/Meta_GWAS/2023_gwas/rg2022/rg22_ch_GCST90165267_buildGRCh38.p5e8.tsv",
                  header=T)

ch_regen$CHROM <- paste0("chr",ch_regen$chr)
ch_regen.df <- ch_regen[,c(23, 3, 3, 1)]
names(ch_regen.df) <- c("chrom", "start", "end","varID")
gr.ch_regen <- as(ch_regen.df, "GenomicRanges")

ch_hg38tohg19 <- liftOver(gr.ch_regen, 
                              hg38_to_19_chain)
ch_hg38tohg19 <- as.data.frame(ch_hg38tohg19)
ch_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",
                              replacement = "", 
                              x = as.character(ch_hg38tohg19$seqnames))
ch_hg38tohg19$POS_hg19 <- ch_hg38tohg19$start
ch_hg38tohg19 <- merge(ch_hg38tohg19[,c(8,9,10)], 
                           ch_regen, 
                           by.x="varID",
                           by.y="name")

ch_hg38tohg19$varID_hg19 <- paste(ch_hg38tohg19$chr_hg19, 
                                  ch_hg38tohg19$POS_hg19, 
                                  ch_hg38tohg19$other_allele, 
                                  ch_hg38tohg19$effect_allele, 
                                  sep=":")
fwrite(ch_hg38tohg19,"/ch_progression/aric/gwas/rg2022/rg22_ch_GCST90165267_buildGRCh19.p5e8.csv", 
       row.names = F, col.names = T, sep = ",")

  # DNMT3A
dnmt_regen <- fread("/Meta_GWAS/2023_gwas/rg2022/rg22_dnmt3a.GCST90165271_buildGRCh38.p5e8.tsv",
                    header=T)

dnmt_regen$CHROM <- paste0("chr",dnmt_regen$chr)
dnmt_regen.df <- dnmt_regen[,c(23, 3, 3, 1)]
names(dnmt_regen.df) <- c("chrom", "start", "end","varID")
gr.dnmt_regen <- as(dnmt_regen.df, "GenomicRanges")

dnmt3a_hg38tohg19 <- liftOver(gr.dnmt_regen, hg38_to_19_chain)
dnmt3a_hg38tohg19 <- as.data.frame(dnmt3a_hg38tohg19)
dnmt3a_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",replacement = "", 
                                   x = as.character(dnmt3a_hg38tohg19$seqnames))
dnmt3a_hg38tohg19$POS_hg19 <- dnmt3a_hg38tohg19$start
dnmt3a_hg38tohg19 <- merge(dnmt3a_hg38tohg19[,c(8,9,10)], 
                           dnmt_regen, 
                           by.x="varID",
                           by.y="name")

dnmt3a_hg38tohg19$varID_hg19 <- paste(dnmt3a_hg38tohg19$chr_hg19, 
                                      dnmt3a_hg38tohg19$POS_hg19, 
                                      dnmt3a_hg38tohg19$other_allele, 
                                      dnmt3a_hg38tohg19$effect_allele, 
                                      sep=":")

fwrite(dnmt3a_hg38tohg19,"/ch_progression/aric/gwas/rg2022/rg22_DNMT3A_GCST90165271_buildGRCh19.p5e8.csv", 
       row.names = F, col.names = T, sep = ",")
  # TET2
tet_regen <- fread("/Meta_GWAS/2023_gwas/rg2022/rg22_tet2_GCST90165281_buildGRCh38.p5e8.tsv",
                   header=T)

tet_regen$CHROM <- paste0("chr",tet_regen$chr)
tet_regen.df <- tet_regen[,c(23, 3, 3, 1)]
names(tet_regen.df) <- c("chrom", "start", "end","varID")
gr.tet_regen <- as(tet_regen.df, "GenomicRanges")

tet_hg38tohg19 <- liftOver(gr.tet_regen, hg38_to_19_chain)
tet_hg38tohg19 <- as.data.frame(tet_hg38tohg19)
tet_hg38tohg19$chr_hg19 <- gsub(pattern = "chr",
                                replacement = "", 
                                x = as.character(tet_hg38tohg19$seqnames))
tet_hg38tohg19$POS_hg19 <- tet_hg38tohg19$start
tet_hg38tohg19 <- merge(tet_hg38tohg19[,c(8,9,10)], 
                           tet_regen, 
                           by.x="varID",
                           by.y="name")

tet_hg38tohg19$varID_hg19 <- paste(tet_hg38tohg19$chr_hg19, 
                                   tet_hg38tohg19$POS_hg19, 
                                   tet_hg38tohg19$other_allele, 
                                   tet_hg38tohg19$effect_allele, 
                                  sep=":")
fwrite(tet_hg38tohg19,"/ch_progression/aric/gwas/rg2022/rg22_TET_GCST90165281_buildGRCh19.p5e8.csv", 
       row.names = F, col.names = T, sep = ",")

### Extract ARIC GWAS
regen.ch_dnmt_tet <- as.data.frame(rbind(ch_hg38tohg19,
                                         dnmt3a_hg38tohg19,
                                         tet_hg38tohg19))
library(dplyr)
regen.ch_dnmt_tet.v2 <- regen.ch_dnmt_tet[, c(1,2,3,8,
                                              11,14,15,
                                              24,26)] %>% 
  group_by(varID_hg19) %>% 
  summarise(varID_hg38=unique(varID),
            chrom=as.numeric(unique(chr_hg19)),
            POS_hg19=unique(POS_hg19),
            Traits = paste(trait, collapse = "|"),
            OR_UKB=max(odds_ratio),
            Pval_UKB=min(p_value),
            AAF_UKB=max(effect_allele_frequency),
            INFO_UKB=paste(additional_info, collapse = "|"))

# fwrite(regen.ch_dnmt_tet.v2,"/ch_progression/aric/gwas/rg2022/rg22_CH_DNMT3A_TET_GCST_buildGRCh19.p5e8.csv", 
#        row.names = F, col.names = T, sep = ",")

  ## Extract ARIC SNPs
setwd("/ch_progression/aric/gwas/gwas_current/meta_aric/")
# while read pheno; do zgrep -wf <(cut -d',' -f1 /ch_progression/aric/gwas/rg2022/rg22_CH_DNMT3A_TET_GCST_buildGRCh19.p5e8.csv | awk 'NR>1' | sort -V | uniq) /ch_progression/aric/gwas/gwas_current/meta_aric/meta_aric_N3015.chr1_22.ea2378_aa637.incident_${pheno}.out.gz > /ch_progression/aric/gwas/gwas_current/meta_aric/regen_GWAS_Cat_variants.${pheno}.tsv; done < <(echo -e "CH\nCH_or_growingClones\nDTA\nSF\nDDR\nDNMT3A\nTET2\nASXL1" ) &
ch.aric <- fread("regen_GWAS_Cat_variants.CH.tsv", header = F, 
                 col.names = c("varID_hg19","EA",	"RefA", "EAF",	
                               "beta",	"se",	"beta_95L",	
                               "beta_95U", "z",	"P",	
                               "log10_P", "q_statistic","q_p",	
                               "i2",	"n_studies",	
                               "n_samples", "effects"))
ch.aric <- merge(ch.aric, regen.ch_dnmt_tet.v2, by="varID_hg19")
qqman::manhattan(ch.aric, chr = "chrom", bp = "POS_hg19", p = "P",snp = "varID_hg19", suggestiveline = -log10(0.05))
plot(-log10(ch.aric$P), -log10(ch.aric$Pval_UKB) )

table(ch.aric$chrom[ch.aric$P<0.05 & ch.aric$EAF>=0.01 & ch.aric$EAF<=0.99])

dnmt3a.aric <- fread("regen_GWAS_Cat_variants.DNMT3A.tsv", header = F, 
                 col.names = c("varID_hg19","EA",	"RefA", "EAF",	
                               "beta",	"se",	"beta_95L",	
                               "beta_95U", "z",	"P",	
                               "log10_P", "q_statistic","q_p",	
                               "i2",	"n_studies",	
                               "n_samples", "effects"))
dnmt3a.aric <- merge(dnmt3a.aric, 
                    regen.ch_dnmt_tet.v2, 
                    by="varID_hg19")
table(dnmt3a.aric$chrom[dnmt3a.aric$P<0.05 & 
                          dnmt3a.aric$EAF>=0.01 & 
                          dnmt3a.aric$EAF<=0.99])

tet.aric <- fread("regen_GWAS_Cat_variants.TET2.tsv", header = F, 
                     col.names = c("varID_hg19","EA",	"RefA", "EAF",	
                                   "beta",	"se",	"beta_95L",	
                                   "beta_95U", "z",	"P",	
                                   "log10_P", "q_statistic","q_p",	
                                   "i2",	"n_studies",	
                                   "n_samples", "effects"))
tet.aric <- merge(tet.aric, 
                  regen.ch_dnmt_tet.v2, 
                  by="varID_hg19")
table(tet.aric$chrom[tet.aric$P<0.05 & 
                    tet.aric$EAF>=0.01 & 
                    tet.aric$EAF<=0.99])

(tet.aric[tet.aric$P<0.05 & 
                 tet.aric$EAF>=0.01 & 
                 tet.aric$EAF<=0.99 & 
                 tet.aric$chrom==14, 
         c("varID_hg19", "P")])
# TET: 14:96188847:T:A P=0.023625
# CH: 14:96175006:A:C P=0.016226
(dnmt3a.aric[dnmt3a.aric$P<0.05 & 
               dnmt3a.aric$EAF>=0.01 & 
               dnmt3a.aric$EAF<=0.99 & 
               dnmt3a.aric$chrom==14, 
          c("varID_hg19", "P")])
########################

dnmt3a.aric.v2 <- merge(dnmt3a.aric, var2.rg22_ch_dnmt_tet, 
                     by ="varID_hg19", all.x=T)

tet.aric.v2 <- merge(tet.aric, var2.rg22_ch_dnmt_tet, 
                        by ="varID_hg19", all.x=T)

ch.aric.v2 <- merge(ch.aric, var2.rg22_ch_dnmt_tet, 
                        by ="varID_hg19", all.x=T)
# Lead variants from RG22 GWAS
require(openxlsx)
list_of_datasets <- list("CH"=ch.aric.v2,"DNMT3A"=dnmt3a.aric.v2,
                         "TET2"=tet.aric.v2)
write.xlsx(list_of_datasets, file = "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableX.rg22_sig_vars.xlsx")

fwrite(dnmt3a.aric.v2, "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/Table.DNMT3A.rg22_sig_vars.csv",
       row.names = F, col.names = T, sep=",", na = "NA")

fwrite(ch.aric.v2, "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/Table.CH.rg22_sig_vars.csv",
       row.names = F, col.names = T, sep=",", na = "NA")
fwrite(tet.aric.v2, "/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/Table.TET2.rg22_sig_vars.csv",
       row.names = F, col.names = T, sep=",", na = "NA")
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==10,c(1,4,5,9)])

(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==3,c(1,4,5,9)])

(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==6,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==7,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==8,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==12,c(1,4,5,9)])

(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==11,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==13,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==14,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==17,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==18,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==20,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==21,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==22,c(1,4,5,9)])

(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==1,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==2,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==4,c(1,4,5,9)])
(var2.rg22_ch_dnmt_tet[var2.rg22_ch_dnmt_tet$chrom==5,c(1,4,5,9)])
####### 
## read modified files
library(readxl)
ch_overall <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableX.rg22_sig_vars.xlsx", 
                         sheet = 2)
sort(table(ch_overall$Loci))

dnmt3a_only <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableX.rg22_sig_vars.xlsx", 
                          sheet = 3)
sort(table(dnmt3a_only$Loci))
tet2_only <- read_excel("/Baylor_ARIC_Exomes/Paper1_ARIC/GWAS/TableX.rg22_sig_vars.xlsx", 
                        sheet = 4)

sort(table(tet2_only$Loci[tet2_only$P<=0.05]))
sort(table(dnmt3a_only$Loci[dnmt3a_only$P<=0.05]))
sort(table(ch_overall$Loci[ch_overall$P<=0.05]))

(table(tet2_only[tet2_only$P<=0.05, c("Loci", "effects")]))

##### forest plot for GWAS
library(data.table)
library(meta)
library(grid)
library(scales)
###
plot_gwas <- fread("/Baylor_ARIC_Exomes/Paper1_ARIC/Display/gwas_table_for_plot_GWAS.p05.csv",
                   header = T, nrows = 28)
plot_gwas$P_val <- formatC(x = plot_gwas$P, digits = 1,format = "E")
plot_gwas$OR <- formatC(round(exp(plot_gwas$beta),2))
plot_gwas$CI95 <- paste0("[",
                         formatC(round(exp(plot_gwas$beta_95L),2), 
                                   digits = 2, format = "f"),
                       ", ",
                       formatC(round(exp(plot_gwas$beta_95U),2), digits = 2, format = "f"), 
                       "]")
plot_gwas$or_95ci <- paste0(plot_gwas$OR," ",plot_gwas$CI95)
plot_gwas$Inicident_Outcome <- ordered(plot_gwas$Inicident_Outcome, 
                                       levels = c("Overall CH", 
                                                  "DNMT3A", "TET2"))
plot_gwas$Loci = ordered(plot_gwas$Loci, 
                         levels=c("SMC4", "TERT", 
                                  "ZNF318", "CD164", "HBS1L",
                                  "H2AFV","GSDMC","ATM",
                                  "TCL1A", "DNAH2", "SETBP1",
                                  "BCL2L1", "RUNX1"))

plot_gwas$var_genes <- paste0("chr",plot_gwas$varID_hg38,
                              " (", plot_gwas$Loci,")")

plot_gwas$var_genes = ordered(plot_gwas$var_genes,
                           levels=c("chr3:160279759:C:T (SMC4)", "chr3:160453665:C:T (SMC4)", "chr3:160541081:A:T (SMC4)", "chr5:1277462:G:A (TERT)", "chr5:1282299:G:A (TERT)", "chr5:1285859:C:A (TERT)", "chr5:1390234:T:C (TERT)", "chr6:109314783:G:A (CD164)", "chr6:109440234:A:G (CD164)", "chr6:135114363:T:C (HBS1L)", "chr6:43428587:C:T (ZNF318)", "chr7:44887228:G:A (H2AFV)", "chr8:129518150:G:C (GSDMC)","chr8:129605779:C:T (GSDMC)", "chr11:108283370:G:C (ATM)", "chr14:95708669:A:C (TCL1A)", "chr14:95722510:T:A (TCL1A)", "chr17:7827723:C:G (DNAH2)", "chr18:44491062:A:G (SETBP1)", "chr18:44558619:A:G (SETBP1)", "chr18:44617198:C:T (SETBP1)", "chr20:31681257:G:A (BCL2L1)", "chr21:34972056:C:T (RUNX1)")) 
plot_gwas$pval_prevoutcome <- paste0(plot_gwas$Pval_UKB," (", 
                                     plot_gwas$Prevalent_outcome, 
                                     ")")
plot_gwas$Empty_space_L <- "     "
plot_gwas$Empty_space_R <- "     "
## Format metagen
b_gwas <- metagen(TE = beta,
                    lower = beta_95L, 
                    upper = beta_95U,
                    studlab = Inicident_Outcome,
                    subgroup=var_genes,
                    data=plot_gwas,
                    sm="OR")


pdf("/Baylor_ARIC_Exomes/Paper1_ARIC/Display/GWAS.Forest_incidentCH.v1.pdf", 
    width = 12, height= 16)
forest(x = b_gwas, 
       sortvar = plot_gwas$var_genes,
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 2), 
       smlab="Effect of SNVs\non Incident CH\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab", "EA", "EAF"), 
       leftlabs = c("                     ","EA","EAF","              "),
       rightcols=c( "OR","CI95","P_val","effects", "Prevalent_outcome"),
       rightlabs=c("OR","95% CI","P","Direction","Prevalent CH (Kessler et al. 2022)"),
       #rightcols=NULL, 
       #rightlabs=NULL,
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)
dev.off()

## v2
b_gwas.v2 <- metagen(TE = beta,
                     lower = beta_95L, 
                     upper = beta_95U,
                     studlab = var_genes,
                     subgroup=Inicident_Outcome,
                     data=plot_gwas,
                     sm="OR")

pdf("/Baylor_ARIC_Exomes/Paper1_ARIC/Display/GWAS.Forest_incidentCH.v2.pdf", 
    width = 13, height= 8.5)
forest(x = b_gwas.v2, 
       sortvar = plot_gwas$var_genes,
       common=F, 
       random=F, 
       hetstat=F, 
       subgroup=k.w>=1, 
       weight.study="same",  
       level=0.95, 
       xlim=c(0.5, 2), 
       smlab="Effect of SNVs\non Incident CH\n", 
       smlab.pos=0, 
       colgap=unit(7, "mm"),
       xlab="Odds Ratio", 
       squaresize=0.6, 
       col.subgroup="black", 
       colgap.left=unit(0.1,"cm"),
       colgap.forest.left="3mm", 
       colgap.forest.right="2mm", 
       leftcols=c("studlab", "EA", "EAF","Empty_space_L"), 
       leftlabs = c("                     ","EA","EAF",""),
       rightcols=c( "Empty_space_R", "or_95ci", "P_val",
                    "effects", "Prevalent_outcome"),
       rightlabs=c("", "OR [95% CI]",
                   "P",   "Direction",
                   "     Prevalent CH\n(Kessler et al. 2022)"),
       just.addcols.right = c("center","center", 
                              "center", "center","left"),
       #rightcols=NULL, 
       #rightlabs=NULL,
       col.inside="black", 
       plotwidth=unit(6.5, "cm"), 
       print.subgroup.name=F)
dev.off()

#############################
########
# hg19 combo with hg38
load(file = "/CHIP_annotation/2022_CHIP_Call/ARIC_hg38/March30.all_hg38_novaseq1_10.rda")
names(pass_AN.chip.hg38_aric) <- make.names(names(pass_AN.chip.hg38_aric),unique = TRUE)

hg19_hg38.chip <- merge(hg19.chip.v2, 
                        pass_AN.chip.hg38_aric, 
                        by="GWAS_ID_Visit")

hg19_hg38.chip_all <- merge(hg19.chip.v2, 
                            pass_AN.chip.hg38_aric, 
                            by="GWAS_ID_Visit", all=T)




########
## Genetic determinant of Clonal Hematopoiesis in ARIC
setwd("/Baylor_ARIC_Exomes/Paper1_ARIC/")

library(data.table)

## egtlgen 
# eqtls <- fread("/ch_progression/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
# il_genes <- eqtls$GeneSymbol[grepl(pattern = "il6|il1|aim2",x = eqtls$GeneSymbol,ignore.case = T)]
# sort(table(il_genes))
  ## Sample size >=20k
# ils <- subset(eqtls, eqtls$GeneSymbol %in% 
   #             c("AIM2", "CASP1","CASP5", 
    #              "IFNGR2", "IL10", "IL18BP",
     #             "IL18RAP", "IL1B", "IL1R1",
      #            "IL1R2", "IL6", "IL6ST",
       #           "JAK2", "NEK7", "NRLC4",
        #          "NLRP3", "TNF", "TYK2",
         #         "CARD8", "IFNGR1", "IL18",
          #        "IL18R1", "IL1RAP", "JAK3", 
           #       "STAT4", "STAT6", "IL6R",   
            #      "IL10RA", "IL10RB") & 
             #   eqtls$NrSamples>=20000)

# zgrep -E 'rs6542082|rs6542081|rs10172772|rs7589655|rs10265117|rs80004026|rs11766947|rs6976090|rs35116860|rs4845626|rs11265618|rs12048091|rs10908839|rs2229238|rs7514452|rs11999802|rs7033052|rs4372063|rs4587378|rs6476941|rs2149556|rs7036034|rs76576833' /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/EurAm.chr{2,7,1,9}_incident_TET2.regenie.gz|cut -f1-6,9-13

###
  # CH phenotypes for GWAS
setwd("/ch_progression/aric/gwas/rg2022/")
# aric_baseline_n_v05 <- fread("aric_baseline_longitudinal.N4189_ch_status.2023Feb20.csv", header = T)
aric_baseline_n_v05 <- fread("aric_baseline_longitudinal.N4189_ch_status_trajectories.2023Feb21.csv", header = T)

  # dbGaP linker
# aric_dbgap_linker <- fread("/dbGap/ARIC/87809/RootStudyConsentSet_phs000280.RootStudy.v7.p1.c1.HMB-IRB/PhenotypeFiles/phs000280.v7.pht001441.v7.p1.ARIC_Sample.MULTI.txt.gz", skip = 10, fill = TRUE, header = T)
# aric_dbgap_linker <- aric_dbgap_linker[,c(1,2,4,5)]
# aric_dbgap_linker <- subset(aric_dbgap_linker, !duplicated(aric_dbgap_linker$SAMPLE_ID))
# table(aric_baseline_n_v05$gwasid %in% aric_dbgap_linker$SAMPLE_ID, exclude = NULL)
# TRUE 
# 4189

## link dbGaP ids with GWASID
# aric_baseline_n_v05 <- merge(aric_dbgap_linker, 
#                              aric_baseline_n_v05,
#                              by.x = "SAMPLE_ID",
#                              by.y = "gwasid")

### PCA:
# plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 --extract <(zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --pca --out ARIC_EurAm_chr1_22.aaray_snps
# plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --extract <(zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --pca --out /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.aray_snps
pca_ea <- fread("/ch_progression/aric/pca/ARIC_EurAm_chr1_22.aray_snps.eigenvec", header = T)
# pca_ea$IID <- as.character(pca_ea$IID)
table(aric_baseline_n_v05$SUBJECT_ID %in% pca_ea$IID)
# FALSE  TRUE 
# 1497  2692
aric_baseline_n_v05_ea <- merge(aric_baseline_n_v05, pca_ea, 
                                by.x="SUBJECT_ID", by.y="IID")

aric_baseline_n_v05_ea$age_base2 <- aric_baseline_n_v05_ea$age_base^2

aric_baseline_n_v05_ea_gwas <- aric_baseline_n_v05_ea[, c(1,1,52,53,2:5,13,128,114,37,38,39,40,10,7,58:84,115,116,118:127)]
names(aric_baseline_n_v05_ea_gwas) <- c("FID","IID", names(aric_baseline_n_v05_ea_gwas)[3:56])

# pca_aa <- fread("/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.aray_snps.eigenvec", header = T)
pca_aa <- fread("/ch_progression/aric/pca/ARIC_AfrAm_chr1_22.aray_snps.eigenvec", header = T)

# pca_aa$IID <- as.character(pca_aa$IID)
# plot(pca_ea$PC1, pca_ea$PC2)
table(aric_baseline_n_v05$SUBJECT_ID %in% pca_aa$IID)
# FALSE  TRUE 
# 3488   701
aric_baseline_n_v05_aa <- merge(aric_baseline_n_v05, pca_aa, 
                                by.x="SUBJECT_ID", by.y="IID")
aric_baseline_n_v05_aa$age_base2 <- aric_baseline_n_v05_aa$age_base^2

aric_baseline_n_v05_aa_gwas <- aric_baseline_n_v05_aa[, c(1,1,52,53,2:5,13,128,114,37,38,39,40,10,7,58:84,115,116,118:127)]
names(aric_baseline_n_v05_aa_gwas) <- c("FID","IID", names(aric_baseline_n_v05_aa_gwas)[3:56])

# fwrite(aric_baseline_n_v05_ea, "aric_baseline_n_v05_ea_PCA.2023Feb20.csv", 
#        row.names = F, col.names = T, sep=",", na="NA")
# fwrite(aric_baseline_n_v05_ea_gwas, "aric_baseline_n_v05_ea_PCA_GWAS.2023Feb20.tsv",
#        row.names = F, col.names = T, sep="\t", na="NA", quote = F)
# 
# fwrite(aric_baseline_n_v05_aa, "aric_baseline_n_v05_aa_PCA.2023Feb20.csv", 
#        row.names = F, col.names = T, sep=",", na="NA")
# fwrite(aric_baseline_n_v05_aa_gwas, "aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv",
#        row.names = F, col.names = T, sep="\t", na="NA", quote = F)
  
  ## 21 Feb 2023
  ## Updated: trajectories
# fwrite(aric_baseline_n_v05_ea, "aric_baseline_n_v05_ea_PCA.2023Feb21.csv",
#        row.names = F, col.names = T, sep=",", na="NA")
# fwrite(aric_baseline_n_v05_ea_gwas, "aric_baseline_n_v05_ea_PCA_GWAS.2023Feb21.tsv",
#        row.names = F, col.names = T, sep="\t", na="NA", quote = F)
# 
# fwrite(aric_baseline_n_v05_aa, "aric_baseline_n_v05_aa_PCA.2023Feb21.csv",
#        row.names = F, col.names = T, sep=",", na="NA")
# fwrite(aric_baseline_n_v05_aa_gwas, "aric_baseline_n_v05_aa_PCA_GWAS.2023Feb21.tsv",
#        row.names = F, col.names = T, sep="\t", na="NA", quote = F)
### Prev CH GWAS P<5e-8
# hg37
tet2 <- fread("/Meta_GWAS/rerun/CHIP_GWAS/rerun/p5e-8/summary_p5e8.TET2.tsv")
dnmt3a <- fread("/Meta_GWAS/rerun/CHIP_GWAS/rerun/p5e-8/summary_p5e8.DNMT3A.tsv")
ch <- fread("/Meta_GWAS/rerun/CHIP_GWAS/rerun/p5e-8/summary_p5e8.CHIP.tsv")

# library(devtools)
# install_github("cgrace1978/manhplot", dependencies = T, force = T)
library(manhplot)
infile<-system.file("extdata","cad.add.160614_manhformat.txt.gz",package = "manhplot")
configfile<-system.file("extdata","config.txt", package = "manhplot")
snpfile<-system.file("extdata","56cad.add.160614.variants.txt", package = "manhplot")

## Run manhattan++ with the default paramaters and files included in the package
manhplusplot(infile = infile,
             outfile = "test", 
             configfile = configfile, 
             snpfile = snpfile)
             
manhplusplot(infile = "/CHIP_GWAS/2023_metagwas/650k/overlap/summary_p5e8.CHIP.tsv",
             outfile = "ch_test",
             showgenes = T,chrname = "CHR", 
             posname = "POS",pvalname = "P", 
             frqname = "EAF")

## incident CH GWAS 
tet2_gwas <- fread(cmd="gzcat /ch_progression/aric/gwas/EurAm.chr1_22_incident_TET2.regenie.tsv.gz | cut -f1,2,5,6,10-13", 
                   header = T, fill=TRUE)

tet2_gwas_prev <- merge(tet2_gwas, tet2, by.x="Name", by.y = "MarkerID") 

tet2_gwas.eqtl <- merge(tet2_gwas, ils, by.x="Name", by.y = "SNP")

tet2p5e4 <- subset(tet2_gwas, tet2_gwas$P<5e-4)

rm(tet2_gwas)

dnmt3a_gwas <- fread(cmd="gzcat /ch_progression/aric/gwas/EurAm.chr1_22_incident_DNMT3A.regenie.tsv.gz | cut -f1,2,5,6,10-13", 
                                  header = T, fill=TRUE)

dnmt3a_gwas_prev <- merge(dnmt3a_gwas, dnmt3a, by.x="Name", by.y = "MarkerID") 

dnmt3a_gwas.eqtl <- merge(dnmt3a_gwas, ils, by.x="Name", by.y = "SNP")

dnmt3ap5e4 <- subset(dnmt3a_gwas, dnmt3a_gwas$P<5e-4)

rm(dnmt3a_gwas)

asxl1_gwas <- fread(cmd="gzcat /ch_progression/aric/gwas/EurAm.chr1_22_incident_ASXL1.regenie.tsv.gz| cut -f1,2,5,6,10-13", 
                   header = T, fill=TRUE)
# asxl1_gwas_prev <- merge(asxl1_gwas, asxl1, by.x="Name", by.y = "MarkerID") 

asxl1_gwas.eqtl <- merge(asxl1_gwas, ils, by.x="Name", by.y = "SNP")

asxl1ap5e4 <- subset(asxl1_gwas, asxl1_gwas$P<5e-4)

rm(asxl1_gwas)

ch_gwas <- fread(cmd="gzcat /ch_progression/aric/gwas/EurAm.chr1_22_incident_CH.regenie.tsv.gz| cut -f1,2,5,6,10-13", 
                 header = T, fill=TRUE)
ch_gwas_prev <- merge(ch_gwas, ch, by.x="Name", by.y = "MarkerID") 

ch_gwas.eqtl <- merge(ch_gwas, ils, by.x="Name", by.y = "SNP")

chap5e4 <- subset(ch_gwas, ch_gwas$P<5e-4)

rm(dnmt3a_gwas)

ch_gwas_prev <- subset(ch_gwas_prev, ch_gwas_prev$P.x<0.05)

barplot(table(ch_gwas_prev$CHR))
plot(y=-log10(ch_gwas_prev$P.x), x=(ch_gwas_prev$CHR))
 

dnmt_gwas_prev <- subset(dnmt3a_gwas_prev, dnmt3a_gwas_prev$P.x<0.05)
barplot(table(dnmt_gwas_prev$CHR))
plot(y=-log10(dnmt_gwas_prev$P.x), x=(dnmt_gwas_prev$CHR))

tet2_gwas_prev <- subset(tet2_gwas_prev, tet2_gwas_prev$P.x<0.05)
barplot(table(tet2_gwas_prev$CHR))
plot(y=-log10(tet2_gwas_prev$P.x), x=(tet2_gwas_prev$CHR))

cat(dnmt3a_gwas.eqtl[(dnmt3a_gwas.eqtl$Name=="rs6542082"),])

### Coloc

# library(coloc)
# tet2_gwas <- fread(cmd="gzcat /ch_progression/aric/gwas/EurAm.chr1_22_incident_TET2.regenie.tsv.gz | cut -f1,2,5,6,13,14,23,24,27", 
  #                 header = T, fill=TRUE)

# tet2_gwas.eqtl <- merge(tet2_gwas, ils[ils$GeneSymbol=="IL1B",], 
  #                      by.x="Name", by.y = "SNP")
 
# res_coloc <- vector("list", nrow(dnmt3a_cpg_list))
# for(i in 1:nrow(dnmt3a_cpg_list)){
# prop.table(table(aric_baseline_n_v05$incident_TET2))
# 
# 0          1 
# 0.95584276 0.04415724 
# qtl_inc_tet2 <- list(
  #  pvalues=tet2_gwas.eqtl$P,
  #  N=tet2_gwas.eqtl$N ,
  #  beta=tet2_gwas.eqtl$BETA,
  #  varbeta=tet2_gwas.eqtl$SE,
  #  snp=tet2_gwas.eqtl$Name,
  #  type="cc", s=0.04415724,
  #  MAF=tet2_gwas.eqtl$AAF)
  
# eqtl.il1 <- list(
 # pvalues=tet2_gwas.eqtl$Pvalue,
 # N=tet2_gwas.eqtl$NrSamples ,
 # zscore=tet2_gwas.eqtl$Zscore,
 # snp=tet2_gwas.eqtl$Name,
 # type="quant",
 # MAF=tet2_gwas.eqtl$AAF)
  
 # res_coloc <- coloc.abf(qtl_inc_tet2, 
 #                           eqtl.il1, 
 #                        p12 = 1e-6)
  
#  dat1 <- runsusie(qtl_inc_tet2)
  
#  res_2 <- coloc.susie(qtl_inc_tet2, 
 #                        eqtl.il1, 
  #                       p12 = 1e-6)
  # test.dnmt3a$summary
# }
# cat(res_coloc[[i]])
# test.dnmt3a < coloc.abf(dnmt3a_dat1, dnmt3a_dat2, p12 = 1e-6)

# library(httr)
# library(dplyr)
# Res_Coloc <- bind_rows(res_coloc) 
# Res_Coloc$cpg <- dnmt3a_cpg_list$Var1
# data.table::fwrite(Res_Coloc, "/Methylation/DNAm/Final_results/coloc_CAD2017.dnmt3a.cis_mQTL_dat.csv", row.names = F, sep=",")


######## WHI GWAS samples
# topmed_all <- fread("/topmed/freeze10/phased/bgen/out/freeze.10b.chr5.pass_only.sample", skip=2, header = F)
## tar -xvf /medpop/esp2/mesbah/datasets/topmed/whi/90370/topmed-dcc/exchange/phs001237_TOPMed_WGS_WHI/Investigator_Data_Products/CHIP_20200127_phs001237.tar.gz
# ch_var_whi_topmed <- fread("/topmed/whi/90370/topmed-dcc/exchange/phs001237_TOPMed_WGS_WHI/Investigator_Data_Products/projects/topmed/downloaded_data/investigator_data_products/20190701_CHIP/processed_files/20200127_CHIP/CHIP_20200127_phs001237/CHIP_variants_20200127_phs001237.txt", 
  #                 header = T)

# ch_whi_topmed <- fread("/topmed/whi/90370/topmed-dcc/exchange/phs001237_TOPMed_WGS_WHI/Investigator_Data_Products/projects/topmed/downloaded_data/investigator_data_products/20190701_CHIP/processed_files/20200127_CHIP/CHIP_20200127_phs001237/CHIP_calls_20200127_phs001237.txt", 
  #                         header = T)
