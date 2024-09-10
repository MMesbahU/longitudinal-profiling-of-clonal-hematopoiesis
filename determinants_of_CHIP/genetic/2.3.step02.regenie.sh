#!/bin/bash

### !!! Note
#### REF allele is the effect allele in MGBB53k GWAS with --ref-first flag
#### To make ALT allele as EA, omit --ref-first flag
### !!!

##########################################################################
source /scripts/useuse
use Anaconda
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
source activate regenie_env
#########################################################################

## April 4, 2023
### updated phenotype
## glm w/o pred: regenie --step 2 --ignore-pred --pgen /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr5.1kg_imp --bt --htp EA --phenoFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --phenoColList incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 --covarFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --covarColList "age_base,age_base_sqr,PC{1:10}" --catCovarList "Sex,Center,v2_vs_other" --verbose --minCaseCount 5 --bsize 400 --minMAC 20 --range 5:753744-1753744 --out chr5.tert.ea
## Firth w/o pred: regenie --step 2 --ignore-pred --pgen /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr5.1kg_imp --bt --htp EA --phenoFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --phenoColList incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 --covarFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --covarColList "age_base,age_base_sqr,PC{1:10}" --catCovarList "Sex,Center,v2_vs_other" --firth --approx --pThresh 0.01 --verbose --minCaseCount 5 --bsize 400 --minMAC 20 --range 5:753744-1753744 --out chr5.tert.ea.firth 1>>tet.log 2>>test.log &
## with pred: regenie --step 2 --pred /ch_progression/aric/gwas/gwas_current/NULL_MODEL.EA_pred.list --pgen /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr5.1kg_imp --bt --htp EA --phenoFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --phenoColList incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 --covarFile /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv --covarColList "age_base,age_base_sqr,PC{1:10}" --catCovarList "Sex,Center,v2_vs_other" --firth --approx --pThresh 0.01 --verbose --minCaseCount 5 --bsize 400 --minMAC 20 --range 5:753744-1753744 --out chr5.tert.ea.firth_pred 1>>tet.pred.log 2>>test.pred.log &


### AfrAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_AA /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/step02.regenie.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/gwas/aric_baseline_n_v05_aa_PCA_GWAS.2023Apr3.v2.tsv incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base_sqr,PC{1:10} Sex,Center,v2_vs_other 400 /ch_progression/aric/gwas/gwas_current/step2 AfrAm 0.01 /ch_progression/aric/gwas/gwas_current/NULL_MODEL.AA_pred.list 20 4

# EurAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_EA /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/step02.regenie.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/gwas/aric_baseline_n_v05_ea_PCA_GWAS.2023Apr3.v2.tsv incident_CH_or_growingClones,incident_CH,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base_sqr,PC{1:10} Sex,Center,v2_vs_other 400 /ch_progression/aric/gwas/gwas_current/step2 EurAm 0.01 /ch_progression/aric/gwas/gwas_current/NULL_MODEL.EA_pred.list 20 4
#############
## 
# AfrAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_AA /clonal_hematopoiesis/longitudinal_CH/aric/gwas/step02.regenie.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base2,PC{1:10} Sex,center,v2_vs_other 400 /ch_progression/aric/gwas AfrAm 0.01 /ch_progression/aric/gwas/NULL_MODEL.AA_pred.list 20 4

# EurAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_EA /clonal_hematopoiesis/longitudinal_CH/aric/gwas/step02.regenie.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/aric_baseline_n_v05_ea_PCA_GWAS.2023Feb20.tsv incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base2,PC{1:10} Sex,center,v2_vs_other 400 /ch_progression/aric/gwas EurAm 0.01 /ch_progression/aric/gwas/NULL_MODEL.EA_pred.list 20 4

## age,age2, dAge
# AfrAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_AfrAm /clonal_hematopoiesis/longitudinal_CH/aric/gwas/step02.regenie.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv incident_CH_or_growingClones,incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base2,dAge,PC{1:10} Sex,center,v2_vs_other 400 /ch_progression/aric/gwas/dAge AfrAm 0.01 /ch_progression/aric/gwas/dAge/NULL_MODEL.AfrAm_pred.list 20 4

# EurAm: qsub -wd /ch_progression/aric/gwas/tmpdir  -t 1-22 -R y -l h_vmem=10G -l h_rt=40:00:00 -pe smp 4 -binding linear:4 -N step2.aric_imp_EA /clonal_hematopoiesis/longitudinal_CH/aric/gwas/step02.regenie.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /ch_progression/aric/aric_baseline_n_v05_ea_PCA_GWAS.2023Feb21.tsv incident_CH_or_growingClones,incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 age_base,age_base2,dAge,PC{1:10} Sex,center,v2_vs_other 400 /ch_progression/aric/gwas/dAge EurAm 0.01 /ch_progression/aric/gwas/dAge/NULL_MODEL.EurAm_pred.list 20 4

####################################################

########################### Input files
pgen_dir=${1} # /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen 

PHENO_FILE=${2} # /ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv

PHENO_COL_LIST=${3} # incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1

COVAR_COL_LIST=${4} # age_base,age_base2,PC{1:10}

CAT_COVAR_COL_LIST=${5} # Sex,center,v2_vs_other

GENO_BLOCK_SIZE=${6} # 400

OUTPUT_PATH=${7} # /ch_progression/aric/gwas

COHORT=${8} # EurAm or AfrAm

P_THRESH=${9} # 0.01

STEP1_LIST_FILE=${10} # /ch_progression/aric/gwas/NULL_MODEL.AA_pred.list

MIN_MAC=${11} # 20

# interval_file=${12}

cpus=${12}

# FAM=${13}

chr=${SGE_TASK_ID}

pgen=${pgen_dir}/ARIC_${COHORT}_chr${chr}.1kg_imp # /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr22.1kg_imp

# chr=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1}' ${interval_file} )
#CHR_RANGE=$(awk -v var=${SGE_TASK_ID} 'NR==var{print $1":"$2}' ${interval_file} ) # 1:1-62312655
#BGEN=${BGEN_PREFIX}.chr${chr}.bgen
#BGEN_SAMPLE=${BGEN_SAMPLE_PREFIX}.chr${chr}.sample
OUT_PREFIX=${OUTPUT_PATH}/${COHORT}.chr${chr}
#########################################################################
######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########
echo -e "##########################\n"
echo -e "Job : ${COHORT} chr${chr}\n"
echo -e "##########################\n"
########

## Run Regenie
regenie \
	--step 2 \
	--pgen ${pgen} \
	--bt \
	--htp ${COHORT} \
	--phenoFile ${PHENO_FILE} \
	--phenoColList ${PHENO_COL_LIST} \
	--covarFile ${PHENO_FILE} \
	--covarColList ${COVAR_COL_LIST} \
	--catCovarList ${CAT_COVAR_COL_LIST} \
	--bsize ${GENO_BLOCK_SIZE} \
	--firth \
	--approx \
	--gz \
	--write-samples \
	--print-pheno \
	--pThresh ${P_THRESH} \
	--pred ${STEP1_LIST_FILE} \
	--minMAC ${MIN_MAC} \
	--threads ${cpus} \
	--verbose \
	--minCaseCount 5 \
	--out ${OUT_PREFIX}

######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

