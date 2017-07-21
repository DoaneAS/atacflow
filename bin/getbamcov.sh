sbam=$1
outprefix=$2

spack load samtools

samtools index ${sbam}


bamCoverage --bam ${sbam} --binSize 20 --outFileFormat bigwig --smoothLength 120 \
            --normalizeUsingRPKM \
            --maxFragmentLength 150 \
            -o ${outprefix}.sizefactors.bw --centerReads --extendReads --numberOfProcessors 12
