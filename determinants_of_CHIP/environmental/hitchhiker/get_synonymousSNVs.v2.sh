#!/bin/bash

## mpos10_45.dp20_400.vaf35
## after finised: to combine all files: 
# zcat /datasets/CHIP/ARIC/hiseq_vcf/baseline/hiseq_baseline.nrow_1_500.mpos10_45.dp20_400.vaf35.tsv.gz | head -1 | gzip -c > /datasets/CHIP/ARIC/hiseq_vcf/baseline/all_HiSeq_baseline.mpos10_45.dp20_400.vaf35.tsv.gz
# for files in $(ls -l /datasets/CHIP/ARIC/hiseq_vcf/baseline/hiseq_baseline.nrow_*.mpos10_45.dp20_400.vaf35.tsv.gz | awk '{print $NF}'); do zcat ${files} | awk 'NR>1{print $0}' | gzip -c >> /datasets/CHIP/ARIC/hiseq_vcf/baseline/all_HiSeq_baseline.mpos10_45.dp20_400.vaf35.tsv.gz; done & 

## grep A00070 /datasets/CHIP/ARIC/hiseq_vcf/baseline/hiseq_baseline.nrow_1_500.tsv | wc -l

source /broad/software/scripts/useuse

use Bcftools


# for i in `seq 1 500 $(wc -l /datasets/CHIP/ARIC/hiseq_vcf/hiseq_vcf_file.list | awk '{print $1}')`; do qsub -R y -wd /datasets/CHIP/ARIC/hiseq_vcf/tmpdir -N var_per_sam_${i}_$(( ${i} + 499 )) -l h_rt=20:00:00 -l h_vmem=20G -pe smp 1 -binding linear:1 /tools/longitudinal-profiling-of-clonal-hematopoiesis/detect_CH/hitchhiker/get_synonymousSNVs.v2.sh /datasets/CHIP/ARIC/hiseq_vcf/hiseq_vcf_file.list ${i} $(( ${i} + 499 )) /datasets/CHIP/ARIC/hiseq_vcf/baseline mpos10_45.dp20_400.vaf35; done


# bcftools view -H -i 'FILTER=="PASS" && MBQ[0]>=30 && MBQ[1]>=30 && MMQ[0]>=60 && MMQ[1]>=60 && FMT/DP>=25 && FMT/AF<=0.3 && FMT/AD[0:0]>=3 && FMT/AD[0:1]>=3  && FMT/F1R2[0:0]>=1 && FMT/F1R2[0:1]>=1 && FMT/F2R1[0:0]>=1 && FMT/F2R1[0:1]>=1 && FMT/SB[0:0]>=1 && FMT/SB[0:1]>=1 && FMT/SB[0:2]>=1 && FMT/SB[0:3]>=1 && ECNT=1 && MPOS>10 && MPOS<45' annot/annot.A09377-filtered.hg38_multianno.vcf.gz | less

## 
#echo -e "Sample_ID\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT" > /datasets/CHIP/ARIC/hiseq_vcf/hiseq_baseline.2024April29.tsv; for vcf_file in $(ls -lh /datasets/CHIP/ARIC/hiseq_vcf/annot/*.vcf.gz | awk '{print $NF}'); do bcftools view -H -i '(FILTER=="PASS" || FILTER=="germline" || FILTER=="weak_evidence" || FILTER=="multiallelic") && ExonicFunc.refGene="synonymous_SNV" && MBQ[0]>=30 && MBQ[1]>=30 && MMQ[0]>=60 && MMQ[1]>=60 && FMT/DP>=25 && FMT/AF<=0.3 && FMT/AD[0:0]>=3 && FMT/AD[0:1]>=3  && FMT/F1R2[0:0]>=1 && FMT/F1R2[0:1]>=1 && FMT/F2R1[0:0]>=1 && FMT/F2R1[0:1]>=1 && FMT/SB[0:0]>=1 && FMT/SB[0:1]>=1 && FMT/SB[0:2]>=1 && FMT/SB[0:3]>=1 && ECNT=1 && MPOS>10 && MPOS<45' ${vcf_file} | awk -v samplesid=$(bcftools query -l ${vcf_file} ) '{print samplesid"\t"$0}' >> /datasets/CHIP/ARIC/hiseq_vcf/hiseq_baseline.2024April29.tsv; done 

## ls -lh /datasets/CHIP/ARIC/hiseq_vcf/annot/*.vcf.gz | awk '{print $NF}' > /datasets/CHIP/ARIC/hiseq_vcf/hiseq_vcf_file.list
##
vcf_list=${1}

num_start=${2}

num_end=${3}

outdir=${4}

filters=${5} # mpos10_45.dp20_400.vaf35
##
echo -e "Sample_ID\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT" > ${outdir}/hiseq_baseline.nrow_${num_start}_${num_end}.${filters}.tsv
## MPOS>10 && MPOS<45
while read vcf_file
	do 
		bcftools view -H -i '(FILTER=="PASS" || FILTER=="weak_evidence") && INFO/ExonicFunc.refGene="synonymous_SNV" && INFO/avsnp150=="." && INFO/AF=="." && MBQ[0]>=30 && MBQ[1]>=30 && MMQ[0]>=60 && MMQ[1]>=60 && MPOS>10 && MPOS<45 && FMT/DP>=20 && FMT/DP<=400 && FMT/AF<=0.35 && FMT/AD[0:0]>=2 && FMT/AD[0:1]>=2  && FMT/F1R2[0:0]>=1 && FMT/F1R2[0:1]>=1 && FMT/F2R1[0:0]>=1 && FMT/F2R1[0:1]>=1 && FMT/SB[0:0]>=1 && FMT/SB[0:1]>=1 && FMT/SB[0:2]>=1 && FMT/SB[0:3]>=1 && ECNT=1' ${vcf_file} | awk -v samplesid=$(bcftools query -l ${vcf_file} ) '{print samplesid"\t"$0}' >> ${outdir}/hiseq_baseline.nrow_${num_start}_${num_end}.${filters}.tsv

	done < <(awk -v nstart=${num_start} -v nend=${num_end} 'NR>=nstart && NR<=nend{print $0}' ${vcf_list}) 
####
gzip -f ${outdir}/hiseq_baseline.nrow_${num_start}_${num_end}.${filters}.tsv

## to combine all 
# zcat /datasets/CHIP/ARIC/hiseq_vcf/baseline/hiseq_baseline.nrow_1501_2000.tsv.gz | head -1 | gzip -c > /datasets/CHIP/ARIC/hiseq_vcf/baseline/all_HiSeq_baseline.tsv.gz

# for files in $(ls -l /datasets/CHIP/ARIC/hiseq_vcf/baseline/hiseq_baseline.nrow_*.tsv.gz | awk '{print $NF}'); do zcat ${files} | awk 'NR>1{print $0}' | gzip -c >> /datasets/CHIP/ARIC/hiseq_vcf/baseline/all_HiSeq_baseline.tsv.gz; done & 

