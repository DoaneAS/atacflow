#!/bin/bash

DIR=$1
Sample=$(basename $DIR)
##slots=$2

if [ -z  ${NSLOTS+x} ]; then NSLOTS=4; fi


spack load pigz
source activate idp2
cd $DIR

#find *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'


echo "Trimming adapter sequences, with command..."
echo "find *_L00*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'"

find *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'

cat *_R1_001.trim.fastq > ${Sample}.R1.trim.fastq
cat *_R2_001.trim.fastq > ${Sample}.R2.trim.fastq

echo "completed trimming"
pigz -p $NSLOTS -c ${Sample}.R1.trim.fastq > ${Sample}.R1.trim.fastq.gz

pigz -p $NSLOTS -c ${Sample}.R2.trim.fastq > ${Sample}.R2.trim.fastq.gz


