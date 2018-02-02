#!/bin/bash -l
#$ -N ATACathena
#$ -j y
#$ -m e
#$ -cwd
#$ -l athena=true
#$ -M ashley.doane@gmail.com
#$ -l h_rt=8:00:00
#$ -pe smp 8
#$ -l h_vmem=2G
#$ -R y
#$ -o /home/asd2007/joblogs

########### SET DEFAULTS ###########
TRIM="YES" #
NUC=0 # run nucleoatac, extended runtime required
ALN="bwa"
ATHENA=1
GENOME="hg38"
##############################

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -f|--folderpath)
            FOLDERPATH="$2"
            shift # past argument
            ;;
        -g|--genome)
            GENOME="$2"
            shift # past argument
            ;;
        -t|--trim)
            TRIM="$2"
            shift # past argument
            ;;
        -a|--aligner)
            ALN="$2"
            shift
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done
echo FOLDER PATH   = "${FOLDERPATH}"
echo GENOME     = "${GENOME}"
echo TRIM READS    = "${TRIM}"
echo ALINGER    = "${ALN}"
if [[ -n $1 ]]; then
    echo "Last line of FILEPATH specified as non-opt/last argument:"
    tail -1 $1
fi

#set -e

# detect kernel then source spack FILEPATH if it exists

#if [ -f /pbtech_mounts/softlib001/apps/EL6/spack/share/spack/setup-env.sh ] ; then
#    . /pbtech_mounts/softlib001/apps/EL6/spack/share/spack/setup-env.sh
#fi
spack load jdk
spack load pigz
spack load bwa
#spack load bowtie2
##source /softlib/apps/EL6/R-v3.3.0/env
spack load bzip2

source activate idp2athena

#if [ -f /home/asd2007/ATACseq/set-env.sh ] ; then
#     . /home/asd2007/ATACseq/set-env.sh




#FOLDERPATH=$1 #FOLDERPATH to all the Samples
#GENOME=$2
#GENOME="hg38" #Number indicating the reference gtf [1 for Human 2 for Mouse 3 for other]

# Uses job array for each Sample in the folder
FILEPATH=$(ls ${FOLDERPATH} | tail -n +${SGE_TASK_ID}| head -1)

#change directory to node directory
cd $TMPDIR

echo "File : $FILEPATH Read from $FOLDERPATH\n"

#Obtain the name of the Sample
Sample=$(basename "$FILEPATH")
Sample=${Sample%%.*}


rsync -r -v -a  $FOLDERPATH/$FILEPATH/*.gz ./
#rsync -r -v -a -z $FOLDERPATH/$FILEPATH/*.sra ./
rsync -r -v -a $FOLDERPATH/$FILEPATH/* ./
#SRA=$(ls *.sra)



#fastq-dump --gzip --split-3 --clip *.sra


mkdir -p ${Sample}

#pigz -p $NSLOTS -c *.R1.trim.fastq > $TMPDIR/${Sample}.R1.trim.fastq.gz
#pigz -p $NSLOTS -c *.R2.trim.fastq > $TMPDIR/${Sample}.R2.trim.fastq.gz
#rm *trim.fastq

#parallel -j ${NSLOTS} 'fastq-dump --split-FILEPATHs {}' ::: *.sra

echo "Processing  $Sample ..."

#Figuring out the reference genome


if [[ $ATHENA == 1 ]] ; then
    REFDIR="/athena/elementolab/scratch/asd2007/reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
    #PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
    export PATH="/home/asd2007/anaconda2/bin:$PATH"
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.conf"
fi
# get genome args

#rsync -avP ${REF} $TMPDIR/#
#rsync -av ${REFbt2} $TMPDIR/

mkdir -p ${TMPDIR}/${Sample}

echo "ls of pwd is"
ls -lrth


#find *_L00*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel 'trimAdapters.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'

#'cutadapt -a adaptors_to_trim -A adaptors_to_trim -q 20 --minimum-length 5 -o {}_R1_cutadapt.fastq.gz -p {}_R2_cutadapt.fastq.gz {}_R1.fastq.gz {}_R2.fastq.gz &> {}.cutadapt'

#cat *R1.trim.fastq > ${Sample}.R1.trim.fastq
#cat *R2.trim.fastq > ${Sample}.R2.trim.fastq

#if [ -f "${TMPDIR}/${Sample}.R1.trim.fastq.gz" ];
#then
#    TRIM=NO
#else
#    TRIM=YES
#fi


spack load pigz
spack load bedtools2

    echo "Trimming adapter sequences, with command..."
    echo "find *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.athena.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'"
    find *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' |parallel -j ${NSLOTS} 'trimAdapters.athena.py -a {}_R1_001.fastq.gz -b {}_R2_001.fastq.gz'
    cat *_R1_001.trim.fastq > ${Sample}.R1.trim.fastq
    cat *_R2_001.trim.fastq > ${Sample}.R2.trim.fastq
    echo "completed trimming"
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R1.trim.fastq > $TMPDIR/${Sample}.R1.trim.fastq.gz
    pigz -p $NSLOTS -c $TMPDIR/${Sample}.R2.trim.fastq > $TMPDIR/${Sample}.R2.trim.fastq.gz
    rsync -av $TMPDIR/${Sample}.R1.trim.fastq.gz ${FOLDERPATH}/${Sample}
    rsync -av $TMPDIR/${Sample}.R2.trim.fastq.gz ${FOLDERPATH}/${Sample}

