#!/bin/bash

# qsub -R y -wd /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/hiseq_pieup -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 -N get_files_n4187 /medpop/esp2/mesbah/tools/longitudinal-profiling-of-clonal-hematopoiesis/detect_CH/get_gcp_pileup.sh /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/hiseq_sample.n4187.tsv /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/U2AF1_hotspot_CH_hiseq.n4187.tsv /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/hiseq_pieup/U2AF1_pileup

source /broad/software/scripts/useuse

use Google-Cloud-SDK

## 
sample_list=${1} # /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/hiseq_sample.n4187.tsv

outFile=${2} # /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/U2AF1_hotspot_CH_hiseq.n4187.tsv

fileDir=${3} # /medpop/esp2/mesbah/datasets/CHIP/baylor/ARIC_CHIP/Baylor_ARIC_Exomes/mpileups/hiseq_pieup/U2AF1_pileup

### Extract pileup for all 4187 samples with follow-up visit
awk 'NR>1{print $2}' ${sample_list} | gsutil -m cp -I ${fileDir}/

echo -e "CHROM\tPOS\tvarID\tREF\tALT\tINFO\tADF\tADR\tDP\tSample\tCRAM_ID" > ${outFile}

for files in $(ls -lh ${fileDir}/*.txt.gz| awk '{print $NF}')
do 
		zcat ${files} | awk -v my_cramid=$(basename ${files} ".U2AF1pileup.txt.gz") 'NR>1 {print $0"\t"my_cramid}' >> ${outFile}
done

gzip ${outFile}


