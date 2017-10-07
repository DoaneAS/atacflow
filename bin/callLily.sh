


PATH=$PATH:/home/asd2007/Tools/UCSC/wigToBigWig




##cat scripts/runLILY.R | R --slave --args $CL.K27ac $OUTDIR 12500 2500 hg19_refseq.ucsc hg19.fa.fai


##SAMPLE=/athena/elementolab/scratch/asd2007/Tools/HMCAN/hmcan/GCBpooled
SAMPLE="/athena/elementolab/scratch/asd2007/Tools/HMCAN/hmcan/atacout/GCBpooled"
OUTDIR="/athena/elementolab/scratch/asd2007/Tools/LILY/GCBpooled"

UCSC="/athena/elementolab/scratch/asd2007/Tools/LILY/src/rose2/rose2/annotation/hg38_refseq.ucsc"
chrz="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"


cat scripts/runLILY.R | R --slave --args $SAMPLE $OUTDIR 12500 2500 $UCSC $chrz



## NB

SAMPLE="/athena/elementolab/scratch/asd2007/Tools/HMCAN/hmcan/atacout/NBpooled"
OUTDIR="/athena/elementolab/scratch/asd2007/Tools/LILY/NBpooled"


cat scripts/runLILY.R | R --slave --args $SAMPLE $OUTDIR 12500 2500 $UCSC $chrz



## to normalize bw files:
cat scripts/renormalizeWig.debug.R | R --slave --args /athena/elementolab/scratch/asd2007/Tools/HMCAN/hmcan atac
