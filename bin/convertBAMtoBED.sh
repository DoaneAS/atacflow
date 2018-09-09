#!/bin/bash
p1=$1

if [ -z "$p1" ]
then
    echo "This converts BAM files to BED files"
    echo "$(basename $0) <bamFile>"
    echo "First input is the bamFile"
    exit
fi

o1=$(echo ${p1} | sed -r 's/\.bam$/.tagAlign.gz/g')

o2=$(echo ${p1} | sed -r 's/\.bam$/.chipEncode.tagAlign.gz/g')

o3=$(echo ${p1} | sed -r 's/\.bam$/.bedpe.gz/g')

o4=$(echo ${p1} | sed -r 's/\.bam$/.tn5.tagAlign.gz/g')


#conda activate atacFlow

bamToBed -i ${p1} | awk 'BEGIN{OFS="\t"} $6=="+" { $2=$2+4; $3=$3 ; $4="N" ; print $0} $6=="-"{ $2=$2; $3=$3-5; $4="N" ; print $0}' | gzip -nc > ${o4}

#bedtools bamtobed -i ${p1} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${o2}

bedtools bamtobed -bedpe -mate1 -i ${p1} | gzip -nc > ${o3}


#zcat ${o2} | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0 }' | gzip -c > ${o4}

#zcat ${pooled_ta_file} |  awk -f $'\t' 'begin {ofs = fs}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${pooled_prefix}.pooled.tn5.tagalign.gz
