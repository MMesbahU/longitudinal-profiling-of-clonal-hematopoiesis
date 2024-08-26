setwd("/medpop/esp2/mesbah/projects/ch_progression/aric/pheno/")
library(data.table)
library(readxl)

##### Visit 05 Samples with ARIC ID
## WES CRAM V05
aric_v05 <- fread("Visit05_ARIC_GWAS_cramid.csv", header = T)

# novaseq_sampinfo <- read_excel("ARIC_CHIP_demographics_Bx001-049_20220329_ForMesbah_sj_Clean.xlsx")
novaseq_sampinfo <- fread("ARIC_CHIP_demographics_Bx001-049_20220329_ForMesbah_sj_Clean.csv", 
                          header=T, na.strings = "")
aric_v05 <- merge(aric_v05,
                  novaseq_sampinfo[novaseq_sampinfo$VisitID=="V05",
                                   c(6,8,9)], by="GWAS_ID")
str(aric_v05)

## ARICID linker
aric_linker_art <- fread("ARIC_TopMED_ID_linker.tsv", header = T)

aric_v05.mapped1 <- merge(aric_v05, aric_linker_art,
                         by.x="GWAS_ID", by.y="gwasid")

names(aric_v05.mapped1) <- c(names(aric_v05.mapped1)[1:6],"ARIC_ID")


## Baseline
aric_baseline <- fread("aric_baseline.N10881_ch_status.2023Mar20.csv", 
                       header = T)

aric_v05.mapped2 <- merge(aric_v05, aric_baseline[,c(1,3)],
                         by.x="GWAS_ID", by.y="gwasid")

names(aric_v05.mapped2) <- c(names(aric_v05.mapped2)[1:6],"ARIC_ID")

aric_v05.mapped <- as.data.frame(rbind(aric_v05.mapped1,
                                       aric_v05.mapped2))

length(unique(aric_v05.mapped$GWAS_ID))
# 4233

## Remove duplicates and save
aric_v05.mapped <- subset(aric_v05.mapped, !duplicated(aric_v05.mapped$GWASID_Visit))
str(aric_v05.mapped)
# fwrite(aric_v05.mapped, "aric_v05.mapped.csv", row.names = F, col.names = T, sep=",")


cat("CH at baseline\n")
ch_var_in_baseline <- fread("../Returned_CH_call/baseline_CH_variants_in_aric_hiseq_novaseq_samevisit.maxDP_noDup.plusPileup.2023Jan30.csv", header = T)
head(ch_var_in_baseline)

cat("CH at Visit 05\n")
ch_var_in_v05 <- fread("../Returned_CH_call/ch_var_in_v05.plusPileup.2023Jan30.csv", header = T)
head(ch_var_in_v05)

## DP>=20; AD>=3, FR/RR/>=1 
## Special Filter for U2AF1: min AD>=5
    # Baseline samples
ch_var_in_baseline_qcd <- subset(ch_var_in_baseline, 
                                 (ch_var_in_baseline$GWAS_ID!="A04109" & ch_var_in_baseline$DP>=20 & ch_var_in_baseline$VAF>=0.02 & 
                                 ch_var_in_baseline$FR.Alt>=1 & ch_var_in_baseline$RR.Alt>=1) &
                                  (ch_var_in_baseline$AD.Alt>=3 & ch_var_in_baseline$Gene !="U2AF1") |
                                 (ch_var_in_baseline$AD.Alt>=5 & ch_var_in_baseline$Gene=="U2AF1")) 
## # based-on IGV review: exclude mutations in "A04109"; 
# 49 indels in same sample; lots of indels in nearby regions as well

str(ch_var_in_baseline_qcd)

sort(table(ch_var_in_baseline_qcd$Gene))

    # Visit 05 samples

ch_var_in_v05_qcd <- subset(ch_var_in_v05, 
                            (ch_var_in_v05$DP>=20 & ch_var_in_v05$VAF>=0.02 & 
                             ch_var_in_v05$FR.Alt>=1 & ch_var_in_v05$RR.Alt>=1) &
                            (ch_var_in_v05$AD.Alt>=3 & ch_var_in_v05$Gene !="U2AF1") |
                            (ch_var_in_v05$AD.Alt>=5 & ch_var_in_v05$Gene=="U2AF1")) 

str(ch_var_in_v05_qcd)

sort(table(ch_var_in_v05_qcd$Gene))

## Save filtered variant files
# fwrite(ch_var_in_v05_qcd, "ch_var_in_v05_qcd.23Mar2023.csv", row.names = F, col.names = T, sep=",", quote = T)
# fwrite(ch_var_in_baseline_qcd, "ch_var_in_baseline_qcd.23Mar2023.csv", row.names = F, col.names = T, sep=",", quote = T)

## Annotate whole baseline file
    # VAF>=2%
aric_baseline$CH_baseline <- ifelse(aric_baseline$gwasid %in% ch_var_in_baseline_qcd$GWAS_ID, 1, 0)
cat("CH Base VAF>=2%:\n")
table(aric_baseline$CH_baseline,exclude = NULL)
    # VAF>=10%
aric_baseline$CHvaf10_baseline <- ifelse(aric_baseline$CH_baseline==1 & aric_baseline$gwasid %in% 
                                           ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$VAF>=0.10], 1, 0)
cat("VAF>=10%:\n")
table(aric_baseline$CHvaf10_baseline,exclude = NULL)

  # CH status in samples with both v2 and v5 WES
aric_baseline$CH_baseline_v05 <- ifelse( aric_baseline$CH_baseline==1 & (aric_baseline$gwasid %in% aric_v05.mapped$GWAS_ID) &
                                        (aric_baseline$gwasid %in% ch_var_in_v05_qcd$GWAS_ID), 1, 
                                         ifelse(aric_baseline$CH_baseline==0 & 
                                                aric_baseline$gwasid %in% aric_v05.mapped$GWAS_ID, 0, NA))
cat("CH present at both visit:\n")
table(aric_baseline$CH_baseline_v05, exclude = NULL)

aric_baseline$CHvaf10_baseline_v05 <- ifelse(aric_baseline$CHvaf10_baseline==1 & 
                                             (aric_baseline$gwasid %in% aric_v05.mapped$GWAS_ID) & 
                                             (aric_baseline$gwasid %in% ch_var_in_v05_qcd$GWAS_ID[ch_var_in_v05_qcd$VAF>=0.10]), 
                                             1, 
                                             ifelse(aric_baseline$CHvaf10_baseline==0 & 
                                                    (aric_baseline$gwasid %in% aric_v05.mapped$GWAS_ID) &
                                                    !(aric_baseline$gwasid %in% ch_var_in_v05_qcd$GWAS_ID[ch_var_in_v05_qcd$VAF>=0.10]), 
                                                    0, NA))
cat("CHvaf10_baseline_v05:\n")
table(aric_baseline$CHvaf10_baseline_v05, exclude = NULL)

## DTA: DNMT3A, TET2, ASXL1
aric_baseline$CH_DTA <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene %in%
                                                              c("DNMT3A", "TET2", "ASXL1") ], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("DTA:\n")
table(aric_baseline$CH_DTA,exclude = NULL)

## Splicing Factor SF: c("SF3B1", "U2AF1", "SRSF2", "ZRSR2")
aric_baseline$CH_SF <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene %in%
                                                              c("SF3B1", "U2AF1", "SRSF2", "ZRSR2") ], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("SF:\n")
table(aric_baseline$CH_SF,exclude = NULL)

## DDR: TP53, PPM1D 
aric_baseline$CH_DDR <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene %in%
                                                              c("PPM1D", "TP53") ], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("DDR:\n")
table(aric_baseline$CH_DDR,exclude = NULL)

aric_baseline$CH_DNMT3A <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene == "DNMT3A"], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("DNMT3A:\n")
table(aric_baseline$CH_DNMT3A,exclude = NULL)

aric_baseline$CH_TET2  <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene =="TET2"], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("TET2:\n")
table(aric_baseline$CH_TET2,exclude = NULL)

aric_baseline$CH_ASXL1  <- ifelse(aric_baseline$CH_baseline==1 & 
                               aric_baseline$gwasid %in% 
                               ch_var_in_baseline_qcd$GWAS_ID[ch_var_in_baseline_qcd$Gene == "ASXL1"], 
                               1, 
                               ifelse(aric_baseline$CH_baseline==0,
                                      0,NA))
cat("ASXL1:\n")
table(aric_baseline$CH_ASXL1,exclude = NULL)

str(aric_baseline)

## Save full baseline file
# fwrite(aric_baseline, "aric_baseline_N10881.pheno_ch_status.23Mar2023.csv", row.names=F, col.names=T, sep=",")

## Vist 05 samples with gwasid, aricid, age, sex
aric_v05.mapped <- fread("aric_v05.mapped.csv", header = T, sep=",")
str(aric_v05.mapped)

# sample w/o baseline wes, coded as "NA"
# AD.ALT >=3; VAF>=2%; F/R>=1; DP>=20
aric_v05.mapped$CH_baseline <- ifelse( (aric_v05.mapped$GWAS_ID %in% aric_baseline$gwasid) & 
                                         (aric_v05.mapped$GWAS_ID %in% 
                                          ch_var_in_baseline_qcd$GWAS_ID), 
                                      1, ifelse( (aric_v05.mapped$GWAS_ID %in% aric_baseline$gwasid) &
                                                !(aric_v05.mapped$GWAS_ID %in% ch_var_in_baseline_qcd$GWAS_ID), 
                                                0, NA )  )


table(aric_v05.mapped$CH_baseline, exclude = NULL)

# 44 V05 samples w/o baseline WES data
# 4189 V05 samples with baseline WES data

## CH at Visit 05
aric_v05.mapped$CH_v05 <- ifelse( (aric_v05.mapped$GWAS_ID %in% aric_baseline$gwasid) & 
                                         (aric_v05.mapped$GWAS_ID %in% 
                                          ch_var_in_v05_qcd$GWAS_ID), 
                                      1, ifelse( (aric_v05.mapped$GWAS_ID %in% aric_baseline$gwasid) &
                                                !(aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID), 
                                                0, NA )  )


table(aric_v05.mapped$CH_v05, exclude = NULL)

## Incident CH at VAF>=2% at V05; absent at baseline
## Both visit available
## Prevalent (ant missing samples) CH coded as "NA" and will be excluded

aric_v05.mapped$incident_CH <- ifelse( (aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) &  
                                         aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                      [!(ch_var_in_v05_qcd$GWAS_ID %in% ch_var_in_baseline_qcd$GWAS_ID)], 1, 
                                       ifelse( ( (aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) & 
                                                aric_v05.mapped$CH_baseline==0 & aric_v05.mapped$CH_v05==0),0, NA))

table(aric_v05.mapped$incident_CH, exclude = NULL)

## Incident CH at VAF>=5% at V05; absent at baseline (VAF<2%) 
aric_v05.mapped$incident_CHvaf05 <- ifelse( (aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) &  
                                              aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID[ ch_var_in_v05_qcd$VAF>=0.05 & !(ch_var_in_v05_qcd$GWAS_ID %in% ch_var_in_baseline_qcd$GWAS_ID)], 1, 
                                            ifelse( ((aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) & 
                                                     (aric_v05.mapped$CH_baseline==0 & aric_v05.mapped$CH_v05==0)), 0, NA))

table(aric_v05.mapped$incident_CHvaf05, exclude = NULL)

## Incident CH at VAF>=10% at V05; absent at baseline (VAF<2%) 
aric_v05.mapped$incident_CHvaf10 <- ifelse( (aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) &  
                                              aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID[ ch_var_in_v05_qcd$VAF>=0.10 & !(ch_var_in_v05_qcd$GWAS_ID %in% ch_var_in_baseline_qcd$GWAS_ID)], 1, 
                                            ifelse( ((aric_v05.mapped$GWAS_ID  %in% aric_baseline$gwasid) & 
                                                     (aric_v05.mapped$CH_baseline==0 & aric_v05.mapped$CH_v05==0)), 0, NA))

table(aric_v05.mapped$incident_CHvaf10, exclude = NULL)

### Gene-categories
# Groupings
aric_v05.mapped$incident_DTA <- ifelse( aric_v05.mapped$incident_CH==1 & 
                                       aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                       [ ch_var_in_v05_qcd$Gene %in% c("ASXL1","DNMT3A","TET2")], 1, 
                                        ifelse( aric_v05.mapped$CH_baseline==0 & 
                                                 !is.na(aric_v05.mapped$CH_baseline) &
                                                   (aric_v05.mapped$CH_v05 ==0 | aric_v05.mapped$incident_CH==0 |  
                                                    aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                                    [!(ch_var_in_v05_qcd$Gene %in% c("ASXL1","DNMT3A","TET2"))]), 
                                               0, NA))
                                                     

table(aric_v05.mapped$incident_DTA, exclude = NULL)

# Splicing factors: SF3B1, U2AF1, SRSF2, ZRSR2
aric_v05.mapped$incident_SF <- ifelse( aric_v05.mapped$incident_CH==1 & 
                                       aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                       [ ch_var_in_v05_qcd$Gene %in% c("SF3B1", "U2AF1", "SRSF2", "ZRSR2")], 1, 
                                        ifelse( aric_v05.mapped$CH_baseline==0 & 
                                                 !is.na(aric_v05.mapped$CH_baseline) &
                                                   (aric_v05.mapped$CH_v05 ==0 | aric_v05.mapped$incident_CH==0 |  
                                                    aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                                    [!(ch_var_in_v05_qcd$Gene %in% c("SF3B1", "U2AF1", "SRSF2", "ZRSR2"))]), 
                                               0, NA))
                                                     

table(aric_v05.mapped$incident_SF, exclude = NULL)

# DNA damage repair: TP53, PPM1D
aric_v05.mapped$incident_DDR <- ifelse( aric_v05.mapped$incident_CH==1 & 
                                       aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                       [ ch_var_in_v05_qcd$Gene %in% c("PPM1D", "TP53")], 1, 
                                        ifelse( aric_v05.mapped$CH_baseline==0 & 
                                                 !is.na(aric_v05.mapped$CH_baseline) &
                                                   (aric_v05.mapped$CH_v05 ==0 | aric_v05.mapped$incident_CH==0 |  
                                                    aric_v05.mapped$GWAS_ID %in% ch_var_in_v05_qcd$GWAS_ID
                                                    [!(ch_var_in_v05_qcd$Gene %in% c("PPM1D", "TP53"))]), 
                                               0, NA))

table(aric_v05.mapped$incident_DDR, exclude = NULL)

## Save full baseline file
# fwrite(aric_v05.mapped, "aric_v05.mapped_N4233.demo_ch_status.23Mar2023.csv", row.names=F, col.names=T, sep=",")

## Vist 05 samples with gwasid, aricid, age, sex
aric_v05.mapped <- fread("aric_v05.mapped_N4233.demo_ch_status.23Mar2023.csv", header = T, sep=",")
str(aric_v05.mapped)

  # longitudinal visit pheno  
aric_longitudinal <- fread("bcm/aric_longitudinal_vanilla.tsv", header = T, sep="\t")
aric_v05_col <- c("aricid", names(aric_longitudinal)[grep(pattern = "v5", x = names(aric_longitudinal), ignore.case = T)])
aric_visit05 <- aric_longitudinal[, ..aric_v05_col]
str(aric_longitudinal)

aric_v05.mapped_pheno <- merge(aric_v05.mapped, 
                               aric_longitudinal[, c(1, 92:109, 
                                                     115, 121, 127,
                                                     133,138,139)], 
                             by.x="ARIC_ID", by.y="aricid")

str(aric_v05.mapped_pheno)

## Scale continus variables
aric_v05.mapped_pheno$bmi_v5_std <- scale(aric_v05.mapped_pheno$bmi_v5)
aric_v05.mapped_pheno$gluc_v5_std <- scale(aric_v05.mapped_pheno$gluc_v5)
aric_v05.mapped_pheno$chol_v5_std <- scale(aric_v05.mapped_pheno$chol_v5)
aric_v05.mapped_pheno$hdl_v5_std <- scale(aric_v05.mapped_pheno$hdl_v5)
aric_v05.mapped_pheno$ldl_v5_std <- scale(aric_v05.mapped_pheno$ldl_v5)
aric_v05.mapped_pheno$tg_v5_std <- scale(aric_v05.mapped_pheno$tg_v5)
aric_v05.mapped_pheno$non_hdl_v5_std <- scale(aric_v05.mapped_pheno$chol_v5 - aric_v05.mapped_pheno$hdl_v5)

str(aric_v05.mapped_pheno)


plot (aric_v05.mapped_pheno$non_hdl_v5_std, 
     (aric_v05.mapped_pheno$chol_v5_std-aric_v05.mapped_pheno$hdl_v5_std))



## Save demo, pheno, and ch_status
# fwrite(aric_v05.mapped_pheno, 
#       "aric_v05.mapped_N4233.phenoV5_demo_ch_status.23Mar2023.csv", 
#       row.names = F, col.names = T, sep=",", na = "NA")


heatmap(cor(aric_v05.mapped_pheno[,c(40,41,42,43,44,45,46)], use="complete"))

aric_v05.mapped_pheno <- fread("aric_v05.mapped_N4233.phenoV5_demo_ch_status.23Mar2023.csv", header=T)
str(aric_v05.mapped_pheno)

aric_baseline <- fread("aric_baseline_N10881.pheno_ch_status.23Mar2023.csv", header=T)

aric_baseline$CH_baseline =NULL

str(aric_baseline)

table(aric_baseline$gwasid %in% aric_v05.mapped_pheno$GWAS_ID, exclude = NULL)

# Merge Baseline and longitudinal data
aric_baseline_n_v05 <- merge(aric_v05.mapped_pheno, aric_baseline, 
                             by.x="GWAS_ID", by.y="gwasid")
str(aric_baseline_n_v05)


table(aric_baseline_n_v05$chol_med_base, aric_baseline_n_v05$chol_med_v5)

table(aric_baseline_n_v05$cig_base, aric_baseline_n_v05$ever_smoke, exclude = NULL)

aric_baseline_n_v05$nonHDL_base_std <- scale(aric_baseline_n_v05$chol_base - aric_baseline_n_v05$hdl_base)
summary(aric_baseline_n_v05$nonHDL_base_std)

aric_baseline_n_v05$hdl_base_std <- scale(aric_baseline_n_v05$hdl_base)
summary(aric_baseline_n_v05$hdl_base_std)

aric_baseline_n_v05$ldl_base_std <- scale(aric_baseline_n_v05$ldl_base)
summary(aric_baseline_n_v05$ldl_base_std)

aric_baseline_n_v05$chol_base_std <- scale(aric_baseline_n_v05$chol_base)
summary(aric_baseline_n_v05$chol_base_std)

aric_baseline_n_v05$tg_base_std <- scale(aric_baseline_n_v05$tg_base)
summary(aric_baseline_n_v05$tg_base_std)

aric_baseline_n_v05$age_base_sqr <- aric_baseline_n_v05$age_base^2
summary(aric_baseline_n_v05$age_base_sqr)

aric_baseline_n_v05$Center <-  as.factor(ifelse(aric_baseline_n_v05$center=="M", "M", 
                                                ifelse(aric_baseline_n_v05$center=="W", "W", "F_J")))

table(aric_baseline_n_v05$Center ) 

# fwrite(aric_baseline_n_v05, "aric_baseline_n_v05_N10881.pheno_ch_status.23Mar2023.csv", 
  #      row.names = F, col.names = T, sep=",", na = "NA")

