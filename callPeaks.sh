#!/bin/bash

#Sample_TH29_3.narrow.p0.1_peaks.narrowPeak


# ========================
# Create pseudoReplicates
# =======================
FINAL_BEDPE_FILE=$1
SAMPLE=$2
SPEC=$3

chrsz="/athena/elementolab/scratch/asd2007/Reference/hg19.chrom.sizes"
# Get total number of read pairs
nlines=$( zcat ${FINAL_BEDPE_FILE} | wc -l  )
nlines=$(( (nlines + 1) / 2  ))

# Shuffle and split BEDPE file into 2 equal parts

# ========================
# Create pseudoReplicates
# =======================



mkdir -p ${SAMPLE}/pseudo_reps


zcat $FINAL_BEDPE_FILE | shuf | split -d -l $((nlines)) - $SAMPLE/pseudo_reps/$SAMPLE.nodup.

# SYS command. line 303

awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' "$SAMPLE/pseudo_reps/$SAMPLE.nodup.00" | \
    gzip -c > $SAMPLE/pseudo_reps/$SAMPLE.nodup.pr1.tagAlign.gz

rm -f "$TMPDIR/$SAMPLE/pseudo_reps/$SAMPLE.nodup.00"


awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' "$SAMPLE/pseudo_reps/$SAMPLE.nodup.01" | \
		gzip -c > gzip -c > $SAMPLE/pseudo_reps/$SAMPLE.nodup.pr2.tagAlign.gz


rm -f "$TMPDIR/$SAMPLE/pseudo_reps/$SAMPLE.nodup.01"
# SYS command. line 308

POOLED_PREFIX="${SAMPLE}/$SAMPLE.nodup"
POOLED_TA_FILE="${SAMPLE}/${SAMPLE}.nodup.tn5.tagAlign.gz"
#POOLED_TA_FILE="${SAMPLE}/${SAMPLE}.nsorted.fixmate.nodup.noM.black.tn5.tagAlign.gz"
PR_PREFIX="$SAMPLE/pseudo_reps/$SAMPLE.nodup"
PR1_TA_FILE="$SAMPLE/pseudo_reps/$SAMPLE.nodup.pr1.tagAlign.gz"
PR2_TA_FILE="$SAMPLE/pseudo_reps/$SAMPLE.nodup.pr2.tagAlign.gz"


zcat ${POOLED_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${POOLED_PREFIX}.pooled.tn5.tagAlign.gz
zcat ${PR1_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${PR_PREFIX}.pr1.tn5.tagAlign.gz
zcat ${PR2_TA_FILE} |  awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${PR_PREFIX}.pr2.tn5.tagAlign.gz
# Get total number of read pairs

# Shuffle and split BEDPE file into 2 equal parts


macs2 callpeak -t   ${POOLED_PREFIX}.pooled.tn5.tagAlign.gz -f BED -n  ${POOLED_PREFIX} -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1
sort -k 8gr,8gr ${POOLED_PREFIX}_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${POOLED_PREFIX}.tn5.narrowPeak.gz &

macs2 callpeak -t   ${PR_PREFIX}.pr1.tn5.tagAlign.gz -f BED -n  ${PR_PREFIX}.pr1 -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1
sort -k 8gr,8gr ${PR_PREFIX}.pr1_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${PR_PREFIX}.pr1.tn5.narrowPeak.gz &

macs2 callpeak -t ${PR_PREFIX}.pr2.tn5.tagAlign.gz -f BED -n  ${PR_PREFIX}.pr2 -g $SPEC  --nomodel --shift -75 --extsize 150 --keep-dup all --bdg --SPMR --call-summits -p 1e-1
sort -k 8gr,8gr ${PR_PREFIX}.pr2_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc > ${PR_PREFIX}.pr2.tn5.narrowPeak.gz

#prefix=Sample_TH29_3
pval_thresh=0.01




##out
#Narrowpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.narrowPeak.gz
#Broadpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.broadPeak.gz
#Gappedpeak file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.gappedPeak.gz
#Fold-enrichment bigWig file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.fc.signal.bw
#-log10(pvalue) bigWig file ${PEAK_OUTPUT_DIR}/${CHIP_TA_PREFIX}.pval.signal.bw
#


zcat ${POOLED_PREFIX}.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${POOLED_PREFIX}.tn5.pval0.1.500k.narrowPeak.gz

zcat ${PR_PREFIX}.pr1.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pr1.tn5.pval0.1.500k.narrowPeak.gz


zcat ${PR_PREFIX}.pr2.tn5.narrowPeak.gz | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pr2.tn5.pval0.1.500k.narrowPeak.gz

POOLED=${POOLED_PREFIX}.tn5.pval0.1.500k.narrowPeak.gz


#zcat ${SAMPLE}/${SAMPLE}.tn5.pf.narrowPeak.gz  | sort -grk8 | head -n 500000 | gzip -c > ${PR_PREFIX}.pooled.tn5.pval0.1.500k.narrowPeak.gz


intersectBed -wo -a <(zcat -f ${POOLED}) -b <(zcat -f ${PR_PREFIX}.pr1.tn5.pval0.1.500k.narrowPeak.gz)\
    | awk 'BEGIN{FS="\t";OFS="\t"} {s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}'  | cut -f 1-10 | sort | uniq \
    | intersectBed -wo -a stdin -b <(zcat -f ${PR_PREFIX}.pr2.tn5.pval0.1.500k.narrowPeak.gz) \
    | awk 'BEGIN{FS="\t";OFS="\t"} {s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | cut -f 1-10 | sort | uniq | gzip -c > ${PR_PREFIX}.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz

# SYS command. line 109


# SYS command. line 112


# SYS command. line 114


# SYS command. line 116

#TASKTIME=$[$(date +%s)-${STARTTIME}]; echo "Task has finished (${TASKTIME} seconds)."








#
#
#
#
#
#
#
#
#
#
#
#
#
#
# idr --samples /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/rep1/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5.pf.pval0.1.500K.narrowPeak.gz /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/rep2/Sample_GCB_138_FAST.R1.trim.PE2SE.nodup.tn5.pf.pval0.1.500K.narrowPeak.gz --peak-list /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/macs2/pooled_rep/Sample_GCB_138.R1.trim.PE2SE.nodup.tn5_pooled.pf.pval0.1.500K.narrowPeak.gz --input-file-type narrowPeak \
#			--output-file /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.unthresholded-peaks.txt --rank p.value --soft-idr-threshold 0.1 \
#			--plot --use-best-multisummit-IDR
#
## SYS command. line 78
#
# idr_thresh_transformed=$(awk -v p=0.1 'BEGIN{print -log(p)/log(10)}')
#
## SYS command. line 81
#
# awk 'BEGIN{OFS="\t"} $12>='"${idr_thresh_transformed}"' {if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0"}' /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.unthresholded-peaks.txt \
#			| sort | uniq | sort -k7n,7n | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.13-col.bed.gz
#
## SYS command. line 84
#
# zcat /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.13-col.bed.gz | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.narrowPeak.gz
#
## SYS command. line 85
#
# zcat /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.13-col.bed.gz | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.12-col.bed.gz
#
## SYS command. line 87
#
# bedtools intersect -v -a /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.13-col.bed.gz -b /home/asd2007/melnick_bcell_scratch/asd2007/Reference/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.13-col.bed.gz
#
## SYS command. line 88
#
# zcat /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.13-col.bed.gz | awk 'BEGIN{OFS="	"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.narrowPeak.gz
#
## SYS command. line 89
#
# zcat /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.13-col.bed.gz | awk 'BEGIN{OFS="	"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | gzip -c > /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.12-col.bed.gz
#
## SYS command. line 91
#
# gzip -f /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.unthresholded-peaks.txt
#
## SYS command. line 92
#
# rm -f /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.13-col.bed.gz /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/asd2007/Projects/DataSets/atacData/atacEC3986/BDS/GCB/Sample_GCB_138/out/peak/idr/true_reps/rep1-rep2/rep1-rep2.IDR0.1.filt.13-col.bed.gz
#
## SYS command. line 94
#
# TASKTIME=$[$(date +%s)-${STARTTIME}]; echo "Task has finished (${TASKTIME} seconds)."
