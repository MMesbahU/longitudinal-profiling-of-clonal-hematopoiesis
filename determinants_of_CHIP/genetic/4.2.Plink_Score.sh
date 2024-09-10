#!/bin/bash

bfile=${1} # /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3
beta_snp=${2} # prev_CH_GWAS.txt
out_prefix=${3} # prs.prev_CH

##
plink2_linux_x86_64_20200831 \
	--bfile ${bfile} \
	--score ${beta_snp} 1 2 3 header sum \
	--out ${out_prefix}

# CH
plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/ea_aric_rg22ch_beta.21SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.ea_aric_rg22ch_beta.21SNV

plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/aa_aric_rg22ch_beta.21SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.aa_aric_rg22ch_beta.21SNV

# DNMT3A
plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/ea_aric_rg22dnmt3a_beta.22SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.ea_aric_rg22dnmt3a_beta.22SNV

plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/aa_aric_rg22dnmt3a_beta.22SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.aa_aric_rg22dnmt3a_beta.22SNV

# TET2
plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_EurAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/ea_aric_rg22tet2_beta.6SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.ea_aric_rg22tet2_beta.6SNV

plink2_linux_x86_64_20200831 --bfile /dbgaps/aric_phs000090_GENEVA/geneva_qsub/aric_pgen/qcd_hapmap3/ARIC_AfrAm_chr1_22.hapmap3 --score /ch_progression/aric/gwas/PRS/aa_aric_rg22tet2_beta.6SNV.txt 12 15 17 header --out /ch_progression/aric/gwas/PRS/prs.aa_aric_rg22tet2_beta.6SNV

