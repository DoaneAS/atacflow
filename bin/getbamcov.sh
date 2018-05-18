sbam=$1
outprefix=$2

spack load samtools

samtools index ${sbam}
source activate deep

bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
            --normalizeUsing RPKM \
            --maxFragmentLength 150 \
            -o ${outprefix}.sizefactors.bw --centerReads --extendReads --numberOfProcessors 8


bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
            --normalizeUsing BPM \
            --maxFragmentLength 150 \
            -o ${outprefix}.sizefactors.bpm.bw --centerReads --extendReads --numberOfProcessors 8


source deactivate
