#!/bin/bash


## remove duplicate snp id (i.e. multi-allelic snps with same snpid)
# for chr in {1..22}; do awk -F'\t' 'NR>1{print $3}' /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_EurAm_chr${chr}.1kg_imp.pvar | sort -V | uniq -c | awk '$1==1{print $2}' >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/chr1_22.noDup.maf1_info30.EA.snplist; done

# for chr in {1..22}; do awk -F'\t' 'NR>1{print $3}' /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr${chr}.1kg_imp.pvar | sort -V | uniq -c | awk '$1==1{print $2}' >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/chr1_22.noDup.maf1_info30.AA.snplist; done
for chr in {1..22}; do awk -F'\t' 'NR>1{print $3}'  /topmed/aric_phs000090_GENEVA/pgen/ARIC_EurAm_chr${chr}.1kg_imp.pvar | sort -V | uniq -c | awk '$1==1{print $2}' >> /topmed/aric_phs000090_GENEVA/pgen/chr1_22.noDup.maf1_info30.EA.snplist; done &


for chr in {1..22}; do awk -F'\t' 'NR>1{print $3}'  /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.1kg_imp.pvar | sort -V | uniq -c | awk '$1==1{print $2}' >> /topmed/aric_phs000090_GENEVA/pgen/chr1_22.noDup.maf1_info30.AA.snplist; done &


## AA
# 
# make bed
# for chr in {1..22}; do /others/tools/plink2_linux_x86_64_20200831 --pfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr${chr}.1kg_imp --extract /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/chr1_22.noDup.maf1_info30.AA.snplist --make-bed --out /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr${chr}.hapmap3; done 1>/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/log_aa_hap3.txt 2>/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/err_aa_hap3.txt &
for chr in {1..22}; do /others/tools/plink2_linux_x86_64_20200831 --pfile /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.1kg_imp --extract /topmed/aric_phs000090_GENEVA/pgen/chr1_22.noDup.maf1_info30.AA.snplist --make-bed --out /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.hapmap3; done 1>>/topmed/aric_phs000090_GENEVA/pgen/tmpdir/log_aa_hap3.txt 2>>/topmed/aric_phs000090_GENEVA/pgen/tmpdir/err_aa_hap3.txt &

#rm /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/aa_list_beds.txt

# for chr in {2..22}; do echo "/dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr${chr}.hapmap3.bed /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr${chr}.hapmap3.bim /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr${chr}.hapmap3.fam" >> /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/aa_list_beds.txt; done
for chr in {2..22}; do echo "/topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.hapmap3.bed /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.hapmap3.bim /topmed/aric_phs000090_GENEVA/pgen/ARIC_AfrAm_chr${chr}.hapmap3.fam" >> /topmed/aric_phs000090_GENEVA/pgen/aa_list_beds.txt; done

## 
use .plink-1.90b
bedDir=/topmed/aric_phs000090_GENEVA/pgen
# plink --bed ${bedDir}/ARIC_AfrAm_chr1.hapmap3.bed --bim ${bedDir}/ARIC_AfrAm_chr1.hapmap3.bim --fam ${bedDir}/ARIC_AfrAm_chr1.hapmap3.fam --merge-list ${bedDir}/aa_list_beds.txt --make-bed --out ${bedDir}/ARIC_AfrAm_chr1_22.hapmap3
plink \
	--bed ${bedDir}/ARIC_AfrAm_chr1.hapmap3.bed \
	--bim ${bedDir}/ARIC_AfrAm_chr1.hapmap3.bim \
        --fam ${bedDir}/ARIC_AfrAm_chr1.hapmap3.fam \
        --merge-list ${bedDir}/aa_list_beds.txt \
	--make-bed \
	--out ${bedDir}/ARIC_AfrAm_chr1_22.hapmap3

# for chr in {1..22}; do /others/tools/plink2_linux_x86_64_20200831 --pfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/ARIC_AfrAm_chr${chr}.1kg_imp --extract /yruan/tools/hapmap3.snplist --maf 0.05 --geno 0.1 --hwe 1e-10 --mind 0.1 --write-snplist allow-dups --write-samples --no-id-header --out /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/snp_2_keep/ARIC_AfrAm_chr${chr}.qcd; done 1>log_aa.txt 2>err_aa.txt &




## EA
for chr in {1..22}; do /others/tools/plink2_linux_x86_64_20200831 --pfile ${bedDir}/ARIC_EurAm_chr${chr}.1kg_imp --extract ${bedDir}/chr1_22.noDup.maf1_info30.EA.snplist --make-bed --out ${bedDir}/ARIC_EurAm_chr${chr}.hapmap3; done 1>>${bedDir}/tmpdir/log_ea_hap3.txt 2>>${bedDir}/tmpdir/err_ea_hap3.txt &
## 
# rm /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ea_list_beds.txt
for chr in {2..22}; do echo "${bedDir}/ARIC_EurAm_chr${chr}.hapmap3.bed ${bedDir}/ARIC_EurAm_chr${chr}.hapmap3.bim ${bedDir}/ARIC_EurAm_chr${chr}.hapmap3.fam" >> ${bedDir}/ea_list_beds.txt; done
# #plink --bed ${bedDir}/ARIC_EurAm_chr1.hapmap3.bed --bim ${bedDir}/ARIC_EurAm_chr1.hapmap3.bim --fam ${bedDir}/ARIC_EurAm_chr1.hapmap3.fam --merge-list ${bedDir}/ea_list_beds.txt --make-bed --out ${bedDir}/ARIC_EurAm_chr1_22.hapmap3
plink \
	--bed ${bedDir}/ARIC_EurAm_chr1.hapmap3.bed \
	--bim ${bedDir}/ARIC_EurAm_chr1.hapmap3.bim \
	--fam ${bedDir}/ARIC_EurAm_chr1.hapmap3.fam \
	--merge-list ${bedDir}/ea_list_beds.txt \
	--make-bed \
	--out ${bedDir}/ARIC_EurAm_chr1_22.hapmap3


### PCA
/others/tools/plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 --extract <(zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.EA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --pca --out /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.aaray_snps

/others/tools/plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --extract <(zcat /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --pca --out /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.aray_snps


## regenie
# source /software/scripts/useuse
use Anaconda
## to install
# conda create -n regenie_env -c conda-forge -c bioconda regenie
## Update from v3.1g to v3.1.3g
# conda update -n regenie_env -c conda-forge -c bioconda regenie
# conda update -n base conda
source activate regenie_env

regenie --step 1 --loocv --print-prs --bed /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --extract <(zcat /ch_progression/aric/chr1_22.AA.maf1_info30.snp2keep.tsv.gz | awk '$6==2{print $1}') --phenoFile /ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv --covarFile /ch_progression/aric/aric_baseline_n_v05_aa_PCA_GWAS.2023Feb20.tsv --phenoColList incident_CH,incident_CHvaf10,incident_DTA,incident_SF,incident_DDR,incident_DNMT3A,incident_TET2,incident_ASXL1 --covarColList age_base,age_base2,PC{1:10} --catCovarList Sex,race_BW,center,v2_vs_other --bt --bsize 400 --lowmem --out /ch_progression/aric/gwas/NULL_MODEL.AA


## prep summary
while read traits; do echo -e "varID\tName\tChr\tPos_hg19\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tINFO\tMAC\tN" > /ch_progression/aric/gwas/AfrAm.chr1_22_${traits}.regenie.tsv; for chr in {1..22}; do zcat /ch_progression/aric/gwas/AfrAm.chr${chr}_${traits}.regenie.gz | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $2":"$3":"$4":"$5"\t"$0"\t"($14+$18)}' >> /ch_progression/aric/gwas/AfrAm.chr1_22_${traits}.regenie.tsv; done; done < <(echo -e "incident_CH\nincident_CHvaf10\nincident_DTA\nincident_SF\nincident_DDR\nincident_DNMT3A\nincident_TET2\nincident_ASXL1") &

while read traits; do echo -e "varID\tName\tChr\tPos_hg19\tREF\tALT\tTrait\tCohort\tModel\tEffect\tLCI_Effect\tUCI_Effect\tP\tAAF\tNum_Cases\tCases_Ref\tCases_Het\tCases_Alt\tNum_Controls\tControls_Ref\tControls_Het\tControls_Alt\tBETA\tSE\tINFO\tMAC\tN" > /ch_progression/aric/gwas/EurAm.chr1_22_${traits}.regenie.tsv; for chr in {1..22}; do zcat /ch_progression/aric/gwas/EurAm.chr${chr}_${traits}.regenie.gz | sed '1d' | sed -e 's:;:\t:g' -e 's:REGENIE_BETA=::g' -e 's:REGENIE_SE=::g' -e 's:INFO=::g' -e 's:MAC=::g' | awk '{print $2":"$3":"$4":"$5"\t"$0"\t"($14+$18)}' >> /ch_progression/aric/gwas/EurAm.chr1_22_${traits}.regenie.tsv; done; done < <(echo -e "incident_CH\nincident_CHvaf10\nincident_DTA\nincident_SF\nincident_DDR\nincident_DNMT3A\nincident_TET2\nincident_ASXL1") &

## Prep summary for locuszoom
zcat /ch_progression/aric/gwas/AfrAm.chr1_22_incident_CH.regenie.tsv.gz | cut -f3-6,13,14,23,24 | awk '$1~/^[0-9]+/' | bgzip -c > /dbgaps/aric_phs000090_GENEVA/locus.AfrAm.chr1_22_incident_CH.regenie.tsv.gz && tabix -s 1 -b 2 -e 2 /dbgaps/aric_phs000090_GENEVA/locus.AfrAm.chr1_22_incident_CH.regenie.tsv.gz  &

