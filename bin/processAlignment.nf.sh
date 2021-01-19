#!/bin/bash

# read from command line the unfiltered and unsortde bam file

p1=$1
BLACK=$2
MYSLOTS=$3
NSLOTS="4"
##NSLOTS="${SLURM_NPROCS}"
#NSLOTS=$(($SLURM_NPROCS -1))



#if [ -z "${NSLOTS}+x"  ] ; then
#    if [ -z "${SLURM_NPROCS}+x"  ] ; then
#        NSLOTS="$((${MYSLOTS}-1))"
#    else
#        NSLOTS="$((${SLURM_NPROCS}-1))"
#    fi
#fi



#BLACK="/home/asd2007/melnick_bcell_scratch/asd2007/Reference/encodeBlack.bed"
# help
if [ -z "$p1"  ]
then
    echo "This will sort bam, remove mitochondria/duplicates and histogram insert size"
    echo "$(basename $0) <bamFile>"
    echo "<samFile>: Required bam input file"
    exit
fi


## workaround for conda activate issue
source /home/asd2007/miniconda3/etc/profile.d/conda.sh
conda activate atacFlow

##source ~/.spackloads.sh


#export R_JAVA_LD_LIBRARY_PATH=${JAVA_HOME}/jre/lib/amd64/server
#export PATH="/home/asd2007/Tools/bedtools2/bin:$PATH"

export PICARD="/home/asd2007/Tools/picard/build/libs/picard.jar"

export PATH="/home/asd2007/Tools/picard/build/libs:$PATH"

#export PATH=$JAVA_HOME/bin:$PATH

#alias picard="java -Xms500m -Xmx3G -jar $PICARD"

#export PATH="/athena/elementolab/scratch/asd2007/Tools/homer/bin:$PATH"
echo "Sorting..."
out1prefix=$(echo $p1 | sed 's/\.bam$//')
out1="${out1prefix}.sorted.bam"
echo ${out1}
#samtools view -u -q 30 Sample_Ly7_pooled_500k.bam | sambamba sort --memory-limit 16GB --nthreads 6 /dev/stdin --out Sample_Ly7_pooled_500k.sorted.bam
samtools view -u -q 30 "${p1}" | sambamba sort --memory-limit 16GB --nthreads ${NSLOTS} --out ${out1} /dev/stdin
samtools index $out1
    # echo "aligning : $TMPDIR/${Sample}.R1.trim.fq ,  $TMPDIR/${Sample}.R2.trim.fq using bwa-mem.."
    # bwa mem -t ${NSLOTS} -M $TMPDIR/BWAIndex/genome.fa $TMPDIR/${Sample}.R1.trim.fastq $TMPDIR/${Sample}.R2.trim.fastq | samtools view -bS - >  $TMPDIR/${Sample}.bam
#fi



#samtools rmdup
echo "Removing duplicates..."
out2=$(echo ${out1} | sed 's/\.bam$/.nodup.bam/')
echo ${out2}
picard MarkDuplicates VERBOSITY=WARNING \
       INPUT=${out1} OUTPUT=${out2} \
       METRICS_FILE="${out2}.dups.log" \
       REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
#picard MarkDuplicates INPUT=Sample_N1.sorted.bam OUTPUT=Sample_N1.sorted.nodup.bam METRICS_FILE="sample.dups.log" REMOVE_DUPLICATES=true
# index
samtools index $out2


out2m=$(echo $out1 | sed 's/\.bam$/.nodup.noM.temp.bam/')

out3=$(echo $out1 | sed 's/\.bam$/.nodup.noM.bam/')

export CHROMOSOMES=$(samtools view -H $out2 | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e 'VN:' | sed 's/SN://' | xargs echo)

samtools view -@ ${NSLOTS} -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $out2 $CHROMOSOMES > ${out3}

#samtools idxstats $out2 | cut -f 1 | grep -v chrM | xargs samtools view -b $out2 > $out2m


#something odd happening
out2mb=$(echo $out1 | sed 's/\.bam$/.no.black.bam/')
#bedtools subtract -A -a $out2m -b $BLACK > $out2mb
# Remove multimapping and improper reads

#samtools view -@ 6 -F 1804 -f 2 -u ${out2m} > ${out3}
samtools index $out3


out4=$(echo $out1 | sed 's/\.bam$/.nodup.noM.black.bam/')

#rmBlack.sh ${out3} ${BLACK}



spack load bedtools2@2.27
#spack load /qr4zqdd ##samtools

bedtools subtract -A -a ${out3} -b $BLACK > ${out4}

samtools index ${out4}

out5=$(echo $out1 | sed 's/\.bam$/.nsorted.nodup.noM.bam/')


sambamba sort --memory-limit 16GB -n -t ${NSLOTS} --tmpdir="${TMPDIR}" --out ${out1prefix}.nsorted.nodup.noM.bam ${out4}

#rm ${out2m}








## make bam with marked dups and generate PBC file for QC
## BELOW for QC only
echo "QC processing..."
echo "Namesort ..."

sambamba sort -n -u --memory-limit 16GB \
         --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${out1prefix}.nsort.bam ${p1}

#samtools sort -n $p1 -o ${out1prefix}.nsort.bam

#sambamba sort --memory-limit 30GB \
#         --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${TMPDIR}/${Sample}/${Sample}.bam ${TMPDIR}/${Sample}/${Sample}.bam

samtools fixmate -r ${out1prefix}.nsort.bam ${out1prefix}.nsort.fixmate.bam
#samtools view -F 1804 -f 2 -u  ${out1prefix}.nsort.fixmate.bam | samtools sort - > ${out1prefix}.filt.srt.bam

samtools view -F 1804 -f 2 -u  ${out1prefix}.nsort.fixmate.bam | sambamba sort --memory-limit 16GB \
                                                                          --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${out1prefix}.filt.srt.bam /dev/stdin

#alias picard="java -Xms500m -Xmx5G -jar $PICARD"

picard MarkDuplicates VERBOSITY=ERROR QUIET=TRUE \
        INPUT=${out1prefix}.filt.srt.bam OUTPUT=${out1prefix}.dupmark.bam METRICS_FILE=${out1prefix}.dup.qc VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

#samtools sort -n ${out1prefix}.dupmark.bam -o ${out1prefix}.srt.tmp.bam

sambamba sort --n --memory-limit 16GB \
         --nthreads ${NSLOTS} --tmpdir ${TMPDIR} --out ${out1prefix}.srt.tmp.bam  ${out1prefix}.dupmark.bam

dupmark_bam="${out1prefix}.srt.tmp.bam"

PBC_QC="${out1prefix}.pbc.qc"

bedtools bamtobed -bedpe -i $dupmark_bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_QC}

rm $dupmark_bam

rm ${out1prefix}.dupmark.bam



# histogram file
for w in 1000 500
do
    picard CollectInsertSizeMetrics I=$out4 O="${out4}.window${w}.hist_data" H="${out4}.window${w}.hist_graph.pdf" W=${w}
done

conda deactivate
