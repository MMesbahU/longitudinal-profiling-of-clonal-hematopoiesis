#!/bin/bash

## Step 1: download phs00090_GENEVA phg000248 from dbgap
## Step2: tar -xvf phg000248.v2.GENEVA_ARIC_EA.genotype-imputed-data.c1.HMB-IRB.tar &
## Step 3: comvert imputed gpdosage to pgen file
###
## mkdir -p /topmed/aric_phs000090_GENEVA/pgen/tmpdir
## /dbgap/phs000090_GENEVA_aric/pgen
## Run job: 

### Nov 21, 2023
# mkdir -p /topmed/aric_phs000090_GENEVA/pgen/tmpdir
# qsub -t 1-22 -R y -wd /topmed/aric_phs000090_GENEVA/pgen/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=20:00:00 -N gprob2pgen /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/1.1.gprob_2_pgen.sh /dbgap/phs000090_GENEVA_aric /topmed/aric_phs000090_GENEVA/pgen /others/tools/plink2_linux_x86_64_20200831 /ch_progression/aric/pca 2


########## Remove SNPs: REF==ALT
# for chr in {1..22}; do awk '$4==$5' /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr${chr}.1kg_imp.pvar >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.ref_alt_same.snp2remove.tsv; done &

# for chr in {1..22}; do awk '$4==$5' /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr${chr}.1kg_imp.pvar >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.ref_alt_same.snp2remove.tsv; done &

## Good SNPs
# for chr in {1..22}; do zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.metrics.gz | awk '$3>=0.01 && $3<=0.99 && $4>=0.3 && ($6==0 || $6==2) {print $0}' >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv; done &

# gzip /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv

# for chr in {1..22}; do zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.metrics.gz | awk '$3>=0.01 && $3<=0.99 && $4>=0.3 && ($6==0 || $6==2) {print $0}' >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv; done &

# gzip /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv
#########################	
	## SNPs to Keep
# zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '{print $1}' > /dbgaps/aric_phs000090_GENEVA/geneva_qsub/chr1_22.AA.maf1_info30.snplist

# zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv.gz | awk '{print $1}' > /dbgaps/aric_phs000090_GENEVA/geneva_qsub/chr1_22.EA.maf1_info30.snplist

#########

### Feb 17, 2023
# qsub -t 1-22 -R y -wd /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/tmdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=20:00:00 -N gprob2pgen /tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/1.1.gprob_2_pgen.sh /dbgaps/aric_phs000090_GENEVA/geneva_qsub /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /others/tools/plink2_linux_x86_64_20200831 

########
geno_dir=${1} # "/dbgaps/aric_phs000090_GENEVA/geneva_qsub"

plink_out_dir=${2} #"/dbgaps/aric_phs000090_GENEVA/aric_pgen"

plink2=${3} # "/others/tools/plink2_linux_x86_64_20200831"

dir_snps_2_keep=${4} # /ch_progression/aric/pca # chr1_22.AA.maf1_info30.snplist
# phg000248 version v2
num=${5} # 2

chr=${SGE_TASK_ID}

### Format .fam file
# for chr in {1..22}; do awk 'NR>2 {print $1,$2,0,0,$4,$5}' ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.sample > ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam; done

# for chr in {1..22}; do awk 'NR>2 {print $1,$2,0,0,$4,$5}' ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.sample > ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.fam; done


### Convert gprob to pgen
#
# European American
# for chr in {1..22}; do awk 'NR>2 {print $1,$2,0,0,$4,$5}' ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.sample > ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam; done

# for chr in {1..22}

# do 
	#####################
	# Format .fam file ##
	####################
	## EA: phg000248.v2.GENEVA_ARIC_EA.genotype-imputed-data.c1
	zcat ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}_c1.sample.gz  | awk 'NR>2 {print $1,$2,0,0,$4,$5}' > ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam
	## AA: phg000248.v2.GENEVA_ARIC_AA.genotype-imputed-data.c1
	zcat ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}_c1.sample.gz | awk 'NR>2 {print $1,$2,0,0,$4,$5}' > ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.fam
	
	###################
	## gprob to pgen ##
	##################
	## EA
	${plink2} \
		--extract ${dir_snps_2_keep}/chr1_22.EA.maf1_info30.snplist \
		--import-dosage ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}_c1.gprobs.gz \
		format=3 ref-first skip0=1 skip1=1 noheader single-chr=${chr} pos-col-num=3 \
		--psam ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam \
		--make-pgen \
		--memory 15000 \
		--out ${plink_out_dir}/ARIC_EurAm_chr${chr}.1kg_imp

	## AA
	${plink2} \
		--extract ${dir_snps_2_keep}/chr1_22.AA.maf1_info30.snplist \
                --import-dosage ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}_c1.gprobs.gz \
		format=3 ref-first skip0=1 skip1=1 noheader single-chr=${chr} pos-col-num=3 \
                --psam ${geno_dir}/phg000248.v${num}.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.fam \
		--make-pgen \
		--memory 15000 \
		--out ${plink_out_dir}/ARIC_AfrAm_chr${chr}.1kg_imp


# done





