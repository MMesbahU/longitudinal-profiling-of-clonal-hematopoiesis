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

#### 1 Dec 2023:
### AA: qsub -wd /ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N grt_step1_aric_AA /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/2.2.step01.qt.regenie.sh /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr1_22.hapmap3 /ch_progression/aric/gwas/aa.snp2keep.null_model.list /ch_progression/aric/gwas/grt_aa.dp20.1Dec2023.tsv "dVAF_by_dT,logdVAF_by_dT" "age_base,age_base_sqr,PC{1:10},is_notHiSeq,is_notMUTECT" "Sex,Center,v2_vs_other" 1000 /ch_progression/aric/gwas/growth_rate 4 AA
# qsub -wd /ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N grt_step1_aric_EA_strict /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/2.2.step01.qt.regenie.sh /topmed/aric_phs000090_GENEVA/pgen/ARIC_EurAm_chr1_22.hapmap3 /ch_progression/aric/gwas/ea.snp2keep.null_model.list /ch_progression/aric/gwas/grt_ea.dp20allAD5FRRR2.1Dec2023.tsv "dVAF_by_dT,logdVAF_by_dT" "age_base,age_base_sqr,PC{1:10},is_notHiSeq,is_notMUTECT" "Sex,Center,v2_vs_other" 1000 /ch_progression/aric/gwas/growth_rate 4 EA_strict 

# qsub -wd /ch_progression/aric/gwas/tmpdir -R y -l h_vmem=10G -l h_rt=20:00:00 -pe smp 4 -binding linear:4 -N grt_step1_aric_AA_strict /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/2.2.step01.qt.regenie.sh /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr1_22.hapmap3 /ch_progression/aric/gwas/aa.snp2keep.null_model.list /ch_progression/aric/gwas/grt_aa.dp20allAD5FRRR2.1Dec2023.tsv "dVAF_by_dT,logdVAF_by_dT" "age_base,age_base_sqr,PC{1:10},is_notHiSeq,is_notMUTECT" "Sex,Center,v2_vs_other" 1000 /ch_progression/aric/gwas/growth_rate 4 AA_strict

####
## 
PLINK_Prefix=${1} # /gwas/mgbb53k/qced.mgbb53k.autosomes
snps_list=${2} # /aric/gwas/aa.snp2keep.null_model.list
PHENO_FILE=${3} # /Meta_GWAS/rerun/mgbb53k/mgbb53k.imp_new.noRel_sk.21Jul2022.tsv
PHENO_COL_LIST=${4} # hasCHIP,hasDNMT3A,hasTET2,hasASXL1,hasPPM1D,hasTP53,hasSF3B1,hasZNF318,hasSRSF2,hasJAK2,hasZBTB33
COVAR_COL_LIST=${5} # PC{1:10},Age_Genotyping,sqrAge_Genotyping
CAT_COVAR_COL_LIST=${6} # Sex,Ancestry_Self_cat,Batch_CHIP_call
BIN_SIZE=${7} # 1000
OUTPUT_PATH=${8} # /gwas/mgb53k
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
	--qt \
	--apply-rint \
	--loocv \
	--print-prs \
	--bed ${PLINK_Prefix} \
	--extract ${snps_list} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bsize ${BIN_SIZE} \
	--lowmem \
	--minCaseCount 5 \
	--threads ${cpus} \
	--verbose \
	--out ${OUTPUT_PATH}/NULL_MODEL.${pop_prefix}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

