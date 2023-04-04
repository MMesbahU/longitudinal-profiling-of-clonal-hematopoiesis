#!/bin/bash

########## Remove SNPs: REF==ALT
# for chr in {1..22}; do awk '$4==$5' /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr${chr}.1kg_imp.pvar >> /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.ref_alt_same.snp2remove.tsv; done &

# for chr in {1..22}; do awk '$4==$5' /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr${chr}.1kg_imp.pvar >> /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.ref_alt_same.snp2remove.tsv; done &

## Good SNPs
# for chr in {1..22}; do zcat /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.metrics.gz | awk '$3>=0.01 && $3<=0.99 && $4>=0.3 && ($6==0 || $6==2) {print $0}' >> /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv; done &

# gzip /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv

# for chr in {1..22}; do zcat /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.metrics.gz | awk '$3>=0.01 && $3<=0.99 && $4>=0.3 && ($6==0 || $6==2) {print $0}' >> /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv; done &

# gzip /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv
#########################	
	## SNPs to Keep
# zcat /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '{print $1}' > /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/chr1_22.AA.maf1_info30.snplist

# zcat /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv.gz | awk '{print $1}' > /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/chr1_22.EA.maf1_info30.snplist

#########

### Feb 17, 2023
# qsub -t 1-22 -R y -wd /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/tmdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=20:00:00 -N gprob2pgen /medpop/esp2/mesbah/tools/clonal_hematopoiesis/longitudinal_CH/aric/gwas/1.1.gprob_2_pgen.sh /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub /broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen /medpop/esp2/mesbah/others/tools/plink2_linux_x86_64_20200831 

########
geno_dir=${1} # "/broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/geneva_qsub"

plink_out_dir=${2} #"/broad/hptmp/mesbah/dbgaps/aric_phs000090_GENEVA/aric_pgen"

plink2=${3} # "/medpop/esp2/mesbah/others/tools/plink2_linux_x86_64_20200831"

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
	## EA
	awk 'NR>2 {print $1,$2,0,0,$4,$5}' ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.sample > ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam
	## AA
	awk 'NR>2 {print $1,$2,0,0,$4,$5}' ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.sample > ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.fam
	
	###################
	## gprob to pgen ##
	##################
	## EA
	${plink2} \
		--extract ${geno_dir}/chr1_22.EA.maf1_info30.snplist \
		--import-dosage ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.gprobs.gz \
		format=3 ref-first skip0=1 skip1=1 noheader single-chr=${chr} pos-col-num=3 \
		--psam ${geno_dir}/phg000248.v1.GENEVA_ARIC_EA.genotype-imputed-data.c1/ARIC_EurAm_chr${chr}.fam \
		--make-pgen \
		--memory 15000 \
		--out ${plink_out_dir}/ARIC_EurAm_chr${chr}.1kg_imp

	## AA
	${plink2} \
		--extract ${geno_dir}/chr1_22.AA.maf1_info30.snplist \
                --import-dosage ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.gprobs.gz \
		format=3 ref-first skip0=1 skip1=1 noheader single-chr=${chr} pos-col-num=3 \
                --psam ${geno_dir}/phg000248.v1.GENEVA_ARIC_AA.genotype-imputed-data.c1/ARIC_AfrAm_chr${chr}.fam \
		--make-pgen \
		--memory 15000 \
		--out ${plink_out_dir}/ARIC_AfrAm_chr${chr}.1kg_imp


# done





