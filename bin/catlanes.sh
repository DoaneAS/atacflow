#!/bin/bash


FOLDERPATH=$1



FILEPATH=$(ls ${FOLDERPATH})

Sample=$(basename "$FOLDERPATH")
#Sample=${Sample%%.*}



cat "${FOLDERPATH}"/*_R1_001.fastq.gz >"${FOLDERPATH}"/${Sample}.R1.fastq.gz

cat "${FOLDERPATH}"/*_R2_001.fastq.gz >"${FOLDERPATH}"/${Sample}.R2.fastq.gz


