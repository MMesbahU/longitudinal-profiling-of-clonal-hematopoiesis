#!/bin/bash

##########################################################################
source /broad/software/scripts/useuse
use Anaconda
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
source activate regenie_env
#########################################################################

## Apr 4, 2023: Updated phenotype
# Note: exclude BMI_Cat: cut -f1-38 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.tsv > /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv  
## AA: 
# qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_AA /medpop/esp2/mesbah/tools/longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aa.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aric_baseline_n_v05_aa_PCA_GWAS.2023Apr3.v2.tsv "incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base_sqr,PC{1:10}" "Sex,Center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/gwas_current 4 AA

## EA: 
# qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_EA /medpop/esp2/mesbah/tools/longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/ea.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv "incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base_sqr,PC{1:10}" "Sex,Center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/gwas_current 4 EA


########################### Input files
## regenie --step 1 --loocv --print-prs --bed /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --extract <(zcat /medpop/esp2/mesbah/projects/ch_progression/aric/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --phenoFile /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv --covarFile /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv --phenoColList incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 --covarColList age_base,age_base2,PC{1:10} --catCovarList Sex,race_BW,center,v2_vs_other --bt --bsize 400 --lowmem --out /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/NULL_MODEL.AA

## zcat /medpop/esp2/mesbah/projects/ch_progression/aric/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}' > /broad/hptmp/mesbah/ch_progression/aric/gwas/aa.snp2keep.null_model.list
## zcat /medpop/esp2/mesbah/projects/ch_progression/aric/chr1_22.EA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}' > /broad/hptmp/mesbah/ch_progression/aric/gwas/ea.snp2keep.null_model.list

## GWAS stratified by AA and EA 
### AA: qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_AA /medpop/esp2/mesbah/tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aa.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv  "incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base2,PC{1:10}" "Sex,center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas 4 AA
### EA: qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_EA /medpop/esp2/mesbah/tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/ea.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_ea_PCA_GWAS.2023Feb20.tsv "incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base2,PC{1:10}" "Sex,center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas 4 EA

## adjusted for age,age2,dAge
## AfrAm: qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_AfrAm /medpop/esp2/mesbah/tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/aa.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb21.tsv  "incident_CH_or_growingClones,incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base2,dAge,PC{1:10}" "Sex,center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/dAge 4 AfrAm

### EurAm: qsub -wd /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N step1_aric_EurAm /medpop/esp2/mesbah/tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/step01.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/ea.snp2keep.null_model.list /medpop/esp2/mesbah/projects/ch_progression/aric/aric_baseline_n_v05_ea_PCA_GWAS.2023Feb21.tsv "incident_CH_or_growingClones,incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1" "age_base,age_base2,dAge,PC{1:10}" "Sex,center,v2_vs_other" 1000 /medpop/esp2/mesbah/projects/ch_progression/aric/gwas/dAge 4 EurAm

## 
PLINK_Prefix=${1} # /broad/hptmp/mesbah/gwas/mgbb53k/qced.mgbb53k.autosomes
snps_list=${2} # /broad/hptmp/mesbah/ch_progression/aric/gwas/aa.snp2keep.null_model.list
PHENO_FILE=${3} # /medpop/esp2/mesbah/projects/Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL_LIST=${4} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${5} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${6} # Sex,Ancestry_Self_cat,Batch_CHIP_call
BIN_SIZE=${7} # 1000
OUTPUT_PATH=${8} # /broad/hptmp/mesbah/gwas/mgb53k
cpus=${9}
pop_prefix=${10} # AA or EA

#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

## Run Regenie
regenie \
	--step 1 \
	--loocv \
	--print-prs \
	--bed ${PLINK_Prefix} \
	--extract ${snps_list} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bt \
	--bsize ${BIN_SIZE} \
	--lowmem \
	--threads ${cpus} \
	--verbose \
	--minCaseCount 5 \
	--out ${OUTPUT_PATH}/NULL_MODEL.${pop_prefix}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

