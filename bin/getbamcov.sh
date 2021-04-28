#!/bin/bash

sbam=$1
outprefix=$2
MYSLOTS=$3

if [ -z "${NSLOTS}+x"  ] ; then
    NSLOTS="${MYSLOTS}"
fi

##spack load samtools

samtools index ${sbam}

#source ~/.spackloads.sh
##source activate deep


bamCoverage --bam ${sbam} --binSize 1 --outFileFormat bigwig --Offset -1 --maxFragmentLength 150 \
            -o ${outprefix}.ins.bw --numberOfProcessors $MYSLOTS


##source deactivate
