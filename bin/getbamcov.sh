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


bamCoverage --bam ${sbam} --binSize 20 --outFileFormat bigwig --smoothLength 60 \
            --normalizeUsing BPM \
            --maxFragmentLength 150 \
            -o ${outprefix}.bpm.sizefactors.bw --centerReads --extendReads --numberOfProcessors $MYSLOTS


##source deactivate
