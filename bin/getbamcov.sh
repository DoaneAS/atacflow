sbam=$1
outprefix=$2
MYSLOTS=$3

if [ -z "${NSLOTS}+x"  ] ; then
    NSLOTS="${MYSLOTS}"
fi

spack load samtools

samtools index ${sbam}
source activate deep


bamCoverage --bam ${sbam} --binSize 5 --outFileFormat bigwig --smoothLength 60 \
            --normalizeUsing BPM \
            --maxFragmentLength 150 \
            -o ${outprefix}.sizefactors.bpm.bw --centerReads --extendReads --numberOfProcessors "${NSLOTS}"


source deactivate
