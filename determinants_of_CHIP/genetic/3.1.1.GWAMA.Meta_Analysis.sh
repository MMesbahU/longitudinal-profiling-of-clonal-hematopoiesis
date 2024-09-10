#!/bin/bash

source /scripts/useuse

use Tabix

############### Task: Worked
## 1. Run Meta-Analysis
## 2. sorted summary stats MAF 0.1 and N_studies>1
###########################

######### Clock time
echo -e "Job started at: $(date)"
Job_START=$(date +%s)
#########

##############  ARIC AA EA meta-analysis
## gwas summary list
# while read traits; do echo -e "/ch_progression/aric/gwas/gwas_current/step2/chr1_22.EurAm.incident_${traits}.regenie.tsv.gz\n/ch_progression/aric/gwas/gwas_current/step2/chr1_22.AfrAm.incident_${traits}.regenie.tsv.gz" > /ch_progression/aric/gwas/gwas_current/step2/gwama_input.incident_${traits}.ARIC_AA_EA.txt; done < <(echo -e "CH\nDTA\nDNMT3A\nTET2\nASXL1\nSF\nDDR\nCH_or_growingClones")
##### 
## qsub: EA: 2377; AA: 637
# N=3,014
# while read traits; do qsub -R y -wd /ch_progression/aric/gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.incident_${traits}.ARIC_AA_EA /longitudinal-profiling-of-clonal-hematopoiesis/determinants_of_CH/genetic/01.GWAMA.Meta_Analysis.sh /ch_progression/aric/gwas/gwas_current/step2/gwama_input.incident_${traits}.ARIC_AA_EA.txt /ch_progression/aric/gwas/gwas_current/meta_aric/meta_aric_N3015.chr1_22.ea2378_aa637.incident_${traits}; done < <(echo -e "CH\nDTA\nDNMT3A\nTET2\nASXL1\nSF\nDDR\nCH_or_growingClones")

###############################################
##########################
gwas_list=${1}

output_prefix=${2}

# meta_sum=${3} # summary tsv file for plotting
# EUR
# Nme	Chr	Pos_hg19	REF	ALT	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	P	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	BETA	SE	INFO	MAC	N	SNPID
# 1-603505	1	603505	A	C	incident_CH	EurAm	ADD-WGR-FIRTH	0.985226	0.525903	1.84572	0.962936	0.0320381	470	442	28	0	1908	1783	125	0	-0.014884	0.320289	0.501243	152.373169	2378	1:603505:A:C

## zcat chr1_22.AfrAm.incident_CH.regenie.tsv.gz | awk '{print $27"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"$26}'|head -2
## SNPID	REF	ALT	AAF	BETA	SE	P	N
## 1:60553:A:C	A	C	0.0422746	-0.385716	0.677219	0.568977	637

### AA
## while read traits; do zcat /ch_progression/aric/gwas/gwas_current/step2/chr1_22.AfrAm.incident_${traits}.regenie.tsv.gz | awk '{print $27"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"$26}' > /gwas/chr1_22.AfrAm.incident_${traits}.tsv; done < <(echo -e "CH\nDTA\nDNMT3A\nTET2\nASXL1\nSF\nDDR\nCH_or_growingClones") &
## EA 
## while read traits; do zcat /ch_progression/aric/gwas/gwas_current/step2/chr1_22.EurAm.incident_${traits}.regenie.tsv.gz | awk '{print $27"\t"$4"\t"$5"\t"$13"\t"$22"\t"$23"\t"$12"\t"$26}' > /gwas/chr1_22.EurAm.incident_${traits}.tsv; done < <(echo -e "CH\nDTA\nDNMT3A\nTET2\nASXL1\nSF\nDDR\nCH_or_growingClones") &

## ## gwas summary list
# while read traits; do echo -e "/gwas/chr1_22.EurAm.incident_${traits}.tsv\n/broad/hptmp/mesbah/gwas/chr1_22.AfrAm.incident_${traits}.tsv" > /ch_progression/aric/gwas/gwas_current/step2/gwama_input.incident_${traits}.ARIC_AA_EA.txt; done < <(echo -e "CH\nDTA\nDNMT3A\nTET2\nASXL1\nSF\nDDR\nCH_or_growingClones")


## Run GWAMA
GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker SNPID \
	--name_n N \
	--name_ea ALT \
	--name_nea REF \
	--name_eaf AAF \
	--name_beta BETA \
	--name_se SE \
	--indel_alleles \
	-o ${output_prefix}

## compress
bgzip -f ${output_prefix}.out


######### Clock time #########
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time: " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
######### Clock time ########

