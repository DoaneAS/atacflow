#!/bin/bash

rbed=$1
Sample=$2
species=$3
##chrsz=$4

##source activate idp2
##. /home/asd2007/.spackloads.sh

##conda activate atacFlow


macs2 callpeak -t ${rbed} -f BED -n ${Sample}.tag.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

cp ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.narrowPeak

narrowpeak.py ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.tn5.narrowPeak

##macs2 callpeak -t  ${rbed} -f BED -n ${Sample}.tag.broad -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

##/home/asd2007/ATACseq/broadpeak.py  ${Sample}.tag.broad_peaks.broadPeak ${Sample}.tn5.broadPeak 

#macs2 callpeak -t ${rbed} -f BED -n ${Sample}.bedpe -g ${species} --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3

#narrowpeak.py ${Sample}.bedpe_peaks.narrowPeak ${Sample}.tn5.narrowPeak


macs2 callpeak -t  ${rbed} -f BED -n ${Sample}.tag.broad -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

broadpeak.py  ${Sample}.tag.broad_peaks.broadPeak ${Sample}.tn5.broadPeak 

##broadpeak.py  ${Sample}.bedpe.broad_peaks.broadPeak ${Sample}.tn5.broadPeak



##  ## bigBed for narrowpeak file

##cat ${CHRSZ} | grep -P 'chr[\dXY]+[ \t]' > ${BIGBED}.chrsz.tmp
##zcat ${PEAK} | sort -k1,1 -k2,2n > ${BIGBED}.tmp
##bedClip ${BIGBED}.tmp ${BIGBED}.chrsz.tmp ${BIGBED}.tmp2
##
##bedToBigBed -type=bed6+4 -as=narrowPeak.as  ${BIGBED}.tmp2 ${BIGBED}.chrsz.tmp ${BIGBED}
##rm -f ${BIGBED}.tmp ${BIGBED}.tmp2 ${BIGBED}.chrsz.tmp
