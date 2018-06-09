#!/bin/bash -l


########### SETTINGS ###########
ATHENA=1
##############################

while [[ $# -gt 1 ]]
do
    key="$1"

    case $key in
        -s|--sample)
            Sample="$2"
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
        --align)
            BT2ALN=YES
            ;;
        *)
            # unknown option
            ;;
    esac
    shift # past argument or value
done
echo SAMPLE    = "${Sample}"
echo GENOME     = "${GENOME}"
echo TRIM READS    = "${TRIM}"
echo BT2 ALN    = "${BT2ALN}"
if [[ -n $1 ]]; then
    echo "Last line of FILEPATH specified as non-opt/last argument:"
    tail -1 $1
fi
#cd /athena/elementolab/scratch/asd2007/Projects/DataSets/atacData/atacGiorgio

##source activate bds_atac

#FILEPATH=$(ls ${FOLDERPATH} | tail -n +${SGE_TASK_ID}| head -1)


##FILEPATH=$(ls ${FOLDERPATH} | tail -n +11 | head -1)

#change directory to node directory
##cd "$FOLDERPATH"




#Obtain the name of the Sample
#file=$PWD #path to all the Samples

# Uses job array for each Sample in the folder
#ID=$3
#file=$(ls ${path} | tail -n +${ID}| head -1)

#change directory to node directory
#cd $TMPDIR
#cd ${file}
echo "QC for $file"

#Obtain the name of the Sample
##Sample=$(basename "$FILEPATH")
##Sample=${Sample%%.*}


spack load samtools
spack load bedtools2
spack load gcc
spack load jdk
spack load bzip2

export PICARDROOT=/home/asd2007/Tools/picard/build/libs

#cd ${file}

ATHENA=1




if [ $ATHENA == 1 ]
then
    REFDIR="/athena/elementolab/scratch/asd2007/reference"
    ANNOTDIR="/athena/elementolab/scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
    export PATH="/home/asd2007/anaconda2/bin:$PATH"
else
    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/reference"
    PICARDCONF="/home/asd2007/Scripts/picardmetrics.conf"
fi


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
if [[ $GENOME == "hg19" ]] ; then
	  REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/local/share/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa"
elif [[ $GENOME == "hg38" ]] ; then
    echo "genome is ${GENOME}"
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
    PICARDCONF="/athena/elementolab/scratch/asd2007/reference/hg38/picardmetrics.conf"
    #PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
    #ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
    #REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
    #ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
	  REF="/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    RG="hg38"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz=/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes
    seq=/athena/elementolab/scratch/asd2007/reference/hg38/seq
    gensz=hs
    bwt2_idx=/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    REF_FASTA=/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    species_browser=hg38
    # data for ATAQC
    TSS_ENRICH=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_gencode_tss_unique.bed.gz
    DNASE=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
    PROM=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz
    ENH=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz
    REG2MAP=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz
    REG2MAP_BED=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz
    ROADMAP_META=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt
elif [[ $GENOME == "hg38_ENCODE" ]]; then
    DNASE_BED="${ANNOTDIR}/${GENOME}/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
    BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
    #PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
    #ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
    #REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
    #ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
	  REF="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
	  REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/hg38_ENCODE.chrom.sizes"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    RG="hg38"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    chrsz=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/hg38.chrom.sizes
    seq=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/seq
    gensz=hs
    bwt2_idx=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    REF_FASTA=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    species_browser=hg38
    # data for ATAQC
    TSS_ENRICH=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_gencode_tss_unique.bed.gz
    DNASE=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
    PROM=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz
    ENH=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz
    REG2MAP=/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_ENCODE_dnase_avg_fseq_signal_formatted.txt.gz
    REG2MAP_BED=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_celltype_compare_subsample.bed.gz
    ROADMAP_META=/athena/elementolab/scratch/asd2007/reference/hg38_ENCODE/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt
elif [[ $GENOME == "mm10" ]]; then
    genome=mm10
	  gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Mus_UCSC_ref.gtf"
	  REF="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex"
    REFbt2="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/Genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index"
    BLACK="/athena/elementolab/scratch/asd2007/reference/mm10-blacklist.bed"
    RG="mm10"
    SPEC="mm"
    REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Mmusculus/mm10/seq/"
    #chrsz="/athena/elementolab/scratch/asd2007/reference/mm10.genome.chrom.sizes"
    fetchChromSizes mm10 > mm10.chrom.sizes
    chrsz= $PWD/mm10.chrom.sizes
    rsync -av /home/asd2007/Scripts/picardmetrics.Mouse.conf ./
else
    echo "genome is hg19"
	  REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
	  REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
    BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    #fetchChromSizes hg19 > hg19.chrom.sizes
    chrsz="/athena/elementolab/scratch/asd2007/reference/hg19.chrom.sizes"
    cp $chrsz $PWD/hg19.chrom.sizes
    RG="hg19"
    SPEC="hs"
    REFGen="/athena/elementolab/scratch/asd2007/local/share/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa"
fi

#################
################
#### ATACqc ##
#################

##ALIGNED_BAM="${PWD}/${Sample}/${Sample}.bam"
WORKDIR="$PWD"
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample
SAMPLE="$Sample"
PBC="$Sample.pbc.qc"
FINAL_BAM="${Sample}.sorted.nodup.noM.black.bam"
FINAL_BED="${Sample}.nodup.tn5.tagAlign.gz"

#F1="${Sample}.R1.trim.fastq.gz"
#F2="${Sample}.R2.trim.fastq.gz"
ALIGNED_BAM="${Sample}.sorted.bam"
WORKDIR=$PWD
OUTDIR="qc"
OUTPREFIX=$Sample
INPREFIX=$Sample


#cp "${Sample}/${Sample}.sorted.bam" "${Sample}/${Sample}.sort.bam"

export PICARD="/home/asd2007/Tools/picard/build/libs/picard.jar"
export PATH="/home/asd2007/Tools/picard/build/libs:$PATH"
#JAVA_HOME=/home/akv3001/jdk1.8.0_05
alias picard="java -Xms500m -Xmx6G -jar $PICARD"


samtools index "${Sample}.sorted.bam"
samtools index "${FINAL_BAM}"


mkdir -p QCmetrics
mkdir -p QCmetrics/raw


if [ -f ${Sample}.picardcomplexity.qc ]; then
    echo "found picarmetrics file"
else
#    picardmetrics run -f $PICARDCONF -o QCmetrics/raw ${Sample}.sorted.bam
    cp QCmetrics/raw/*.EstimateLibraryComplexity.log QCmetrics/${Sample}.picardcomplexity.qc
    cp QCmetrics/raw/*.EstimateLibraryComplexity.log ${Sample}.picardcomplexity.qc
fi;



PEAKS="${Sample}/peaks/*tag.broad_peaks.broadPeak"

#getFrip.sh "${FINAL_BAM}" "${PEAKS}" 


#python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



export PATH=/athena/elementolab/scratch/asd2007/local/bin:$PATH

#singularity exec /athena/elementolab/scratch/asd2007/projectshg38/datasets/ngi/image/analysis.img python "${DIR}"/run_ataqc.athena.py
python "${DIR}"/run_ataqc.athena.py --workdir $PWD \
    --outdir qc \
    --outprefix ${Sample} \
    --genome ${GENOME} \
    --ref ${REF} --tss $TSS_ENRICH \
    --dnase ${DNASE} \
    --blacklist ${BLACK} \
    --prom $PROM \
    --enh ${ENH} \
   --reg2map ${REG2MAP} \
   --meta ${ROADMAP_META} \
   --alignedbam "${Sample}.sorted.bam"  \
   --alignmentlog "${Sample}.align.log" \
   --coordsortbam "${Sample}.sorted.bam" \
   --duplog "${Sample}.dup.qc" \
   --pbc "${Sample}.pbc.qc" \
   --finalbam "${FINAL_BAM}" \
   --finalbed "${FINAL_BED}" \
   --bigwig "$Sample.sizefactors.bw" \
   --peaks "${Sample}.tn5.broadPeak.gz" \
   --naive_overlap_peaks "${Sample}.tn5.broadPeak.gz" \
   --idr_peaks "${Sample}.tn5.broadPeak.gz"  --processes 4
  # --naive_overlap_peaks "pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz" \
  # --idr_peaks "peaks/${Sample}.tn5.pf.narrowPeak.gz"  --processes 4
#"IDR/${Sample}.IDR.txt.IDR0.1.filt.narrowPeak.gz" \
   






##
##
##
##
##
##
##
##if [ $ATHENA = 1 ]
##then
##    REFDIR="/athena/elementolab/scratch/asd2007/Reference"
##    ANNOTDIR="/athena/elementolab/scratch/asd2007/Reference"
##    PICARDCONF="/home/asd2007/Scripts/picardmetrics.athena.conf"
##else
##    REFDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
##    ANNOTDIR="/zenodotus/dat01/melnick_bcell_scratch/asd2007/Reference"
##fi
##
##if [ $gtf_path = 1 ]
##then
##    #gtf="/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/akv3001/mm10_UCSC_ref.gtf"
##	REF="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex"
##	REFbt2="${REFDIR}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index"
##  # REFbt="/home/asd20i07/dat02/asd2007/Reference/Homo_sapiens/UCSC/mm10/Sequence/BowtieIndex/genome"
##  BLACK="${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
##  chrsz="/athena/elementolab/scratch/asd2007/Reference/hg19.chrom.sizes"
##  cp $chrsz $PWD/hg19.chrom.sizes
##  RG="hg19"
##  SPEC="hs"
##  REFGen="/athena/elementolab/scratch/asd2007/bin/bcbio/genomes/Hsapiens/hg19/seq/"
##fi
##
##export PICARD="/home/asd2007/Tools/picard/build/libs/picard.jar"
##export PATH="/home/asd2007/Tools/picard/build/libs:$PATH"
###JAVA_HOME=/home/akv3001/jdk1.8.0_05
##alias picard="java -Xms500m -Xmx6G -jar $PICARD"
##
##
## #   rm picardmetrics.conf
##
###    picardmetrics run -f "/home/asd2007/picardmetrics.athena.conf"  -o ${Sample}/QCmetrics/filtered ${Sample}/${Sample}.sorted.nodup.noM.bam
###cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
###cp ${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log ${Sample}/${Sample}.picardcomplexity.qc
##
##
##if [ -f ${Sample}/${Sample}.picardcomplexity.qc ]; then
##    echo "found picarmetrics file"
##else
##    picardmetrics run -f $PICARDCONF -o ${Sample}/QCmetrics/raw ${Sample}/${Sample}.sorted.bam
##    cp ${Sample}/QCmetrics/raw/*.EstimateLibraryComplexity.log ${Sample}/QCmetrics/${Sample}.picardcomplexity.qc
##    cp ${Sample}/QCmetrics/raw/*.EstimateLibraryComplexity.log ${Sample}/${Sample}.picardcomplexity.qc
##fi;
##
##
##
##
###samtools flagstat  ${Sample}/${Sample}.sorted.bam > ${Sample}/QCmetrics/raw/${Sample}.sorted.flagstat.txt
###samtools flagstat $TMPDIR/${Sample}/${Sample}.sorted.nodup.noM.bam > $TMPDIR/${Sample}/QCmetrics/filtered/${Sample}.sorted.nodup.noM.flagstat.txt
##
##
###cp $TMPDIR/${Sample}/QCmetrics/filtered/*.EstimateLibraryComplexity.log $TMPDIR/${Sample}/${Sample}.picardcomplexity.qc
##
##
##
##ALIGNED_BAM="${PWD}/${Sample}/${Sample}.bam"
##WORKDIR=$PWD
##OUTDIR="qc"
##OUTPREFIX=$Sample
##INPREFIX=$Sample
##GENOME='hg19' # This is the only genome that currently works
##
##SAMPLE="$Sample"
### Annotation files
##
###ANNOTDIR="/home/asd2007/melnick_bcell_scratch/asd2007/Reference"
##ANNOTDIR="/athena/elementolab/scratch/asd2007/Reference"
##
##DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
##BLACKLIST_BED="${ANNOTDIR}/${GENOME}/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
##TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
##REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
##PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
##ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
##REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
##ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
##PBC="${WORKDIR}/$SAMPLE.pbc.qc"
##
##
##GENOME='hg19' # This is the only genome that currently works
##SAMPLE="$Sample"
##
##PBC="${WORKDIR}/$SAMPLE.pbc.qc"
##FINAL_BAM="${Sample}/${Sample}.sorted.nodup.noM.bam"
##FINAL_BED="${Sample}/${Sample}.nodup.tn5.tagAlign.gz"
##
##F1="${Sample}.R1.trim.fastq.gz"
##F2="${Sample}.R2.trim.fastq.gz"
##ALIGNED_BAM="${Sample}/${Sample}.sorted.bam"
##WORKDIR=$PWD
##OUTDIR="qc"
##OUTPREFIX=$Sample
##INPREFIX=$Sample
##GENOME='hg19' # This is the only genome that currently works
##
### Annotation files
##X="/athena/elementolab/scratch/asd2007/Reference/hg19/"
##DNASE_BED="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
##BLACKLIST_BED="${ANNOTDIR}/${GENOME}/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
##TSS_BED="${ANNOTDIR}/${GENOME}/hg19_RefSeq_stranded.bed.gz"
##REF_FASTA="${ANNOTDIR}/${GENOME}/encodeHg19Male.fa"
##PROM="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
##ENH="${ANNOTDIR}/${GENOME}/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
##REG2MAP="${ANNOTDIR}/${GENOME}/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
##ROADMAP_META="${ANNOTDIR}/${GENOME}/eid_to_mnemonic.txt"
##
##
###cp "${Sample}/${Sample}.sorted.bam" "${Sample}/${Sample}.sort.bam"
##
##
###python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \
##
### condisder changing peaks to peaks/broadpeak
##
##python /home/asd2007/ATACseq/run_ataqc.athena.py --workdir $PWD/${Sample} \
##    --outdir qc \
##    --outprefix ${Sample} \
##    --genome hg19 \
##    --ref "${REFDIR}/hg19/encodeHg19Male/encodeHg19Male.fa" \
##    --tss "${REFDIR}/hg19/hg19_RefSeq_stranded.bed.gz" \
##    --dnase "${REFDIR}/hg19/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz" \
##    --blacklist "${REFDIR}/hg19/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \
##    --prom "${REFDIR}/hg19/reg2map_honeybadger2_dnase_prom_p2.bed.gz" \
##    --enh "${REFDIR}/hg19/reg2map_honeybadger2_dnase_enh_p2.bed.gz" \
##    --reg2map ${REG2MAP} \
##    --meta ${ROADMAP_META} \
##    --fastq1 ${F1} \
##    --fastq2 ${F2} \
##    --alignedbam ${ALIGNED_BAM} \
##    --alignmentlog "${Sample}/${Sample}.align.log" \
##    --coordsortbam "${Sample}/${Sample}.sorted.bam" \
##    --duplog "${Sample}/${Sample}.dup.qc" \
##    --pbc "${Sample}/${Sample}.pbc.qc" \
##    --finalbam "${FINAL_BAM}" \
##    --finalbed "${FINAL_BED}" \
##    --bigwig "$Sample/$Sample.smooth150.center.extend.fpkm.max150.bw" \
##    --peaks "${Sample}/${Sample}.nodup.tn5.pval0.1.500k.narrowPeak.gz" \
##    --naive_overlap_peaks "${Sample}/pseudo_reps/${Sample}.nodup.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz" \
##    --idr_peaks "${Sample}/IDR/${Sample}.IDR.txt.IDR0.1.filt.narrowPeak.gz"
##
