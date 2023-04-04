#!/bin/bash

source /broad/software/scripts/useuse

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

##### 

# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2021/topmed2021_has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.has${traits}.ukb250k.INFO.tsv.gz\n/medpop/esp2/mesbah/datasets/allofus/gwas_ch/${traits}_ch.chr1_22.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2")

# mkdir -p /broad/hptmp/mesbah/ch_gwas/tmpdir

# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/ch_gwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T74UK200k250kAoU98kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k.txt /broad/hptmp/mesbah/ch_gwas/metaGWAS.${traits}.topmed74k_ukb200k250k_aou98k_mgbb53k_biovu54k; done < <(echo -e "CHIP\nDNMT3A\nTET2")


###############################################

## EUR UKB 200k, 250k GWAMA for finemapping
# while read traits; do qsub -R y -wd /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N ukb200k_250k.eur_metaGWAS.${traits} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/ukb200k_ukb250k.meta_${traits}.EUR.txt /medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/eur_only_summary/ukb200k_ukb250k.meta_${traits}.EUR.gwama /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.ukb200k_ukb250k.meta_${traits}.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2")
###
## EUR GWAS
# echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/eur_topmed/topmed2019_EURonly_samples.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.hasCHIP.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.hasCHIP.ukb250k_WB.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k_EUR.hasCHIP.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_hasCHIP_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasCHIP.EUR_only.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/inputUKB/eurUKB/lifted_hg38.chr1_22.has${traits}.21Aug2021.UKB200k_onlyEURsamples.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/ukb450k/250k/hg38/lifted_hg38.chr1_22.has${traits}.ukb250k_WB.INFO.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/MGBB53k/chr1_22.MGBB53k_EUR.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.EUR_only.UKB200kUKB250kMGBB53kBioVU54k.txt; done < <(echo -e "DNMT3A\nTET2")

	# Submit Qsub
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N eur_metaGWAS.${traits}.ukb_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.has${traits}.EUR_only.UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${traits}.ukb_mgb_biovu /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.${traits}.ukb_mgb_biovu.tsv; done < <(echo -e "DNMT3A\nTET2")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N eur_metaGWAS.${traits}.topmed_ukb_mgb_biovu /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.GWAMA.Meta_Analysis.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/input_files/gwama_input.hasCHIP.EUR_only.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_metaGWAS.${traits}.topmed_ukb_mgb_biovu /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/eur_meta_summary.${traits}.topmed_ukb_mgb_biovu.tsv; done < <(echo -e "CHIP")

################### Multi-Ancestry
### Aug 8, 2022
# TOPMed 2019 ukb200k ukb250k mgbb53k biovu54k
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/ukb450k/250k.step2/ukb250k_gwama_input/lifted_hg38.chr1_22.has${traits}.ukb250k.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T64UK200_250kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/meta_summary.${traits}.TOPMed2019UKB200kUKB250kMGBB53kBioVU54k.tsv; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

# TOPMed UKB450k mgbb53k biovu 54k 
#while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/broad/hptmp/mesbah/gwas/ukb450k/450k.step2/lifted_hg38.chr1_22.has${traits}.ukb450k.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB450kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/ukb450k/metagwas/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.T64UK450kMGB53BVU54 /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB450kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/ukb450k/metagwas/results/metaGWAS.${traits}.TOPMed2019UKB450kMGBB53kBioVU54k; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")

## ##################
## 26 July 2022
# while read traits; do echo -e "/medpop/esp/mesbah/GWAS_CHIP/topmed_GWAS/topmed2020/topmed2019.all_samples.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/inputUKB/allUKB/lifted_hg38.chr1_22.has${traits}.11Aug2021_ukb200k_allsamples.INFO.tsv.gz\n/broad/hptmp/mesbah/gwas/mgbb53k/mgbb53k_gwama.chr1_22.has${traits}.tsv.gz\n/medpop/esp/mesbah/GWAS_CHIP/BioVU/BioVU.saige_has${traits}_results_merged_subset.tsv.gz" > /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt; done < <(echo -e "CHIP\nDNMT3A\nTET2\nASXL1")
## Jobs killed with h_vmem=20G
# while read traits; do qsub -R y -wd /broad/hptmp/mesbah/gwas/mgbb53k/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=40G -l h_rt=10:00:00 -N metaGWAS.${traits}.TOPMed2019UKB200kMGBB53kBioVU54k /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/02.GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.Meta-Analysis/01.PrepSummary/gwama_input.has${traits}.TOPMed2019UKB200kMGBB53kBioVU54k.txt /broad/hptmp/mesbah/gwas/mgbb53k/metaGWAS/metaGWAS.${traits}.TOPMed2019UKB200kMGBB53kBioVU54k; done < <(echo -e "CHIP\nDNMT3A\nTET2")
##########################

############# old
## EUR 3: while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_EUR3Samp.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.euronly_UKB_MGBB_BioVU.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.meta_${genename}.euronly_UKB_MGBB_BioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=10:00:00 -N gwama_allEUR.chip /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.EUR_only_CHIP_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.EUR_only_CHIP_GWAS.TopMedUKbbMGbbBioVU
## All: while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_allSamples.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allsamples.topmed_ukb_mgb_biovu.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/meta4.topukbmgbbiovuSaige/GWAMA.allsamples.meta_${genename}.TopMedUKbbMGbbBioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## Mixed: all + EUR
# while read genename; do qsub -R y -wd /broad/hptmp/mesbah/tmpdir -pe smp 1 -binding linear:1 -l h_vmem=20G -l h_rt=10:00:00 -N gwama_mixedSamples.${genename} /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/01.run.GWAMA.sh /medpop/esp2/mesbah/tools/CHIP_metaAnalysis/03.MetaAnalysis/02.Metal_MetaAnalysis/01.run_GWAMA/Input_GWAMA.allTOPMed_eurUKB_MGBB_BioVU.${genename}_GWAS.txt /broad/hptmp/mesbah/ukb_chip/meta_gwas/eurGWAS_meta4/GWAMA.mixedSamples.meta_${genename}.allTOPmed_EUR_UKB_MGB_BioVU; done < <(echo -e "CHIP\nDNMT3A\nTET2")

## 
gwas_list=${1}

output_prefix=${2}

# meta_sum=${3} # summary tsv file for plotting

## Run GWAMA
GWAMA=/medpop/esp2/mesbah/tools/GWAMA_v2.2.2/GWAMA 
${GWAMA} \
	-i ${gwas_list} \
	-qt \
	--name_marker varID \
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

