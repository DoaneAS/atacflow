#!/bin/bash

BAM=$1
BLACK=$2

if [ -z "$BAM"  ]
then
    echo "This will remove blacklisted regions"
    echo "$(basename $0) <bamfile>"
    echo "<bamFile>: Required bam input file"
    exit
fi

out1prefix=$(echo $BAM | sed 's/\.bam$//')
out1="${out1prefix}.black.bam"

if [ -z "$BLACK" ] ; then
    echo "could not find black list regions, using hg38"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
fi


spack load bedtools2
spack load samtools



bedtools subtract -A -a $BAM -b $BLACK > $out1

samtools index $out1


