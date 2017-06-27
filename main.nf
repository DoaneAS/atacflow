#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */ 
       //params.fragmentLength = 200

//params.index = '/home/asd2007/Scripts/nf/fripflow/sindex.tsv'


    params.index = 'indexly7.tsv'
   //params.index = 'indexpooled.tsv'
params.genome = 'hg38'
params.blacklist = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
genome = file(params.genome)
index = file(params.index)
params.chrsz = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"

params.ref = "/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

////// Print parameters ///////
log.info ''
log.info ''
log.info 'A T A C - S e q  ~ P R O C E S S'
log.info '---------------------------------'
log.info ''
log.info "Index File             : ${params.index}"
log.info ''

index = file(params.index)

bwaref = file(params.ref)
results_path = "$PWD/results"

blacklist=file(params.blacklist)


// Clear pipeline.db file

////// Check input parameters //////

if (!params.index) {
  exit 1, "Please specify the input table file"
}


sizefactors = Channel.from(1)



       /*
 *atacs = Channel
 *       .from(index.readLines())
 *       .map { line ->
 *       def list = line.split()
 *              def bed = file(list[0])
 *              def peaks = file(list[1])
 *              def sname = list[2]
 *              def dprefix = file(list[3])
 *       println bed
 *       println peaks
 *       [ sname, bed, peaks, dprefix ]
 *}
 */



/*
 *fastqs = Channel
 *       .from(index.readLines())
 *       .map { line ->
 *              def list = line.split()
 *              def mergeId = list[0]
 *              def id = list[1]
 *              def path = file(list[2])
 *              def message = '[INFO] '
 *              log.info message
 *              def quality = fastq(path).qualityScore()
 *              [ mergeId, id, path, quality ]
 *    }
*
*
*pooled = Channel
*       .from(index.readLines())
*       .map { line ->
*              def list = line.split()
*              def Sample = list[0]
*              def path = file(list[1])
*              def sizefactor = list[2]
*              def message = '[INFO] '
*              log.info message
*              [ Sample, path, sizefactor ]
*    }
*
*/



fastq = Channel
       .from(index.readLines())
       .map { line ->
              def list = line.split()
              def Sample = list[0]
              def path = file(list[1])
              def reads = file("$path/*{1,2}*.fastq.gz")
              // def readsp = "$path/*{R1,R2}.trim.fastq.gz"
              //  def R1 = file(list[2])
              //    def R2 = file(list[3])
              def message = '[INFO] '
              log.info message
              [ Sample, path, reads ]
    }






//    Channel
//       .fromFilePairs( params.reads )                                             
//       .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
//        .set { read_pairs }
  





process bwamem {


    executor 'sge'
    clusterOptions '-l h_vmem=1G -pe smp 4 -l h_rt=1:00:00 -l athena=true'
    //  cpus 12
        scratch "${TMPDIR}"
            //     publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true


        input:
        set Sample, file(path), file(reads) from fastq
        file(bwaref) from bwaref

        output:
        set Sample, file("${Sample}.bam") into newbam

        script:
            //def reads = "$path/*{R1,R2}.trim.fastq.gz"
        """
        set -o pipefail
        bwa mem -t 4 -M /athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta $reads | samtools view -bS - > ${Sample}.bam
        """

    }






process processbam {

    cpus 18
    scratch "${TMPDIR}"
    publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true

    input:
    set Sample, file(nbam) from newbam
    file(BLACK) from blacklist

    output:
    set Sample, file("${Sample}.sorted.bam") into sortedbam
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbam
    file("*.pbc.qc") into pbcqc
    file("*nsort.fixmate.bam") into fixmatebam
    file("*window500.hist_data") into hist_data
    file("*window500.hist_graph.pdf") into fragsizes

    shell:
    '''
    processAlignment.nf.sh !{nbam} !{BLACK} !{task.cpus}
    '''
}







hist_data.subscribe { println "Received: " + file(hist_data)}

fragsizes.subscribe { println "Received: " + file(fragsizes)}


process nsortbam {
    // tag "$Sample"

        cpus 16
        memory '32 GB'

        scratch "${TMPDIR}"
        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true

        input:
        set Sample, file(finalbam) from finalbam

        output:
        set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
        set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
            //val sf into sizefactors

        script:
            // def sf = 1
            //def fbam = file(finalbam)
            // finalbam.renameTo("${Sample}.sorted.nodup.noM.black.bam")

        """
        samtools index ${finalbam}
        sambamba sort --memory-limit 30GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbam}

       ## mv ${finalbam} ${Sample}.sorted.nodup.noM.black.bam
        samtools index ${Sample}.sorted.nodup.noM.black.bam
        """
}





process bam2bed {

    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true

    input:
    set Sample, file(nsbam) from nsortedbam

    output:
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbed
    set Sample, file("${Sample}.nodup.bedpe.gz") into finalbedpe


    script:
    """
    samtools fixmate ${nsbam} ${Sample}.nsorted.fixmate.nodup.noM.bam
    convertBAMtoBED.sh ${Sample}.nsorted.fixmate.nodup.noM.bam
    cp ${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz  ${Sample}.nodup.tn5.tagAlign.gz
    cp ${Sample}.nsorted.fixmate.nodup.noM.bedpe.gz ${Sample}.nodup.bedpe.gz
    """
        }




process callpeaks {

    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true

    input:
    set Sample, file(rbed) from finalbed

    output:
    set Sample, file("${Sample}.tn5.narrowPeak.gz") into narrowpeak
    set Sample, file("${Sample}.tag.narrow_summits.bed") into summits
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeak

    script:
    """
    macs2 callpeak -t ${rbed} -f BED -n ${Sample}.tag.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3
    /home/asd2007/ATACseq/narrowpeak.py ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.tn5.narrowPeak 

    macs2 callpeak -t  ${rbed} -f BED -n ${Sample}.tag.broad -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1

    /home/asd2007/ATACseq/broadpeak.py  ${Sample}.tag.broad_peaks.broadPeak ${Sample}.tn5.broadPeak 
    """
        }


process signalTrack {
    // tag "$Sample"

    cpus 16

        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false

        input:
        set Sample, file(sbam) from bamforsignal
        val sz from sizefactors


        output:
        set Sample, file("${Sample}.sizefactors.bw") into signal

        script:
        """
        samtools index ${sbam} &&
        bamCoverage --bam ${sbam} --binSize 5 \
            --outFileFormat bigwig --smoothLength 100 \
            --normalizeUsingRPKM --scaleFactor ${sz} \
            --maxFragmentLength 150 \
            -o ${Sample}.sizefactors.bw --centerReads --extendReads --numberOfProcessors ${task.cpus}

        """
        }

/*
*process nsortbam {
*    // tag "$Sample"
*
*        cpus 16
*        memory '32 GB'
*
*        scratch "${TMPDIR}"
*        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true
*
*        input:
*        set Sample, file(sbam) from sortedbam
*
*        output:
*        set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
*        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
*        set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
*        val sf into sizefactors
*
*        script:
*            //def fbam = file(finalbam)
*            // finalbam.renameTo("${Sample}.sorted.nodup.noM.black.bam")
*
*        """
*        samtools index ${finalbam}
*        sambamba sort --memory-limit 30GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbam}
*
*        mv ${finalbam} ${Sample}.sorted.nodup.noM.black.bam
*        samtools index ${Sample}.sorted.nodup.noM.black.bam
*        """
*}
*
*
*
  * *process fixmatebam {
  * *
  * *    // publishDir  "$results_path/$sname", mode: 'move', overwrite: false
  * *
  * *        input:
  * *        set Sample, file(nsbam) from nsortedbam
  * *
  * *        output:
  * *        set Sample, file("${Sample}.nsorted.fixmate.nodup.noM.bam") into fixbam
  * *
  * *        script:
  * *        """
  * *        samtools fixmate ${nsbam} ${Sample}.nsorted.fixmate.nodup.noM.bam
  * *         """
  * *        }
  * *
*
*
*process bam2bed {
*
*    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true
*
*    input:
*    set Sample, file(nsbam) from nsortedbam
*
*    output:
*    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbed
*    set Sample, file("${Sample}.nodup.bedpe.gz") into finalbedpe
*
*
*    script:
*    """
*    samtools fixmate ${nsbam} ${Sample}.nsorted.fixmate.nodup.noM.bam
*    convertBAMtoBED.sh ${Sample}.nsorted.fixmate.nodup.noM.bam
*    cp ${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz  ${Sample}.nodup.tn5.tagAlign.gz
*    cp ${Sample}.nsorted.fixmate.nodup.noM.bedpe.gz ${Sample}.nodup.bedpe.gz
*    """
*        }
*
*
*
*
*process callpeaks {
*
*    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true
*
*    input:
*    set Sample, file(rbed) from finalbed
*
*    output:
*    set Sample, file("${Sample}.tn5.narrowPeak.gz") into narrowpeak
*    set Sample, file("${Sample}.tag.narrow_summits.bed") into summits
*    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeak
*
*    script:
*    """
*    macs2 callpeak -t ${rbed} -f BED -n ${Sample}.tag.narrow -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --call-summits -p 1e-3
*    /home/asd2007/ATACseq/narrowpeak.py ${Sample}.tag.narrow_peaks.narrowPeak ${Sample}.tn5.narrowPeak 
*
*    macs2 callpeak -t  ${rbed} -f BED -n ${Sample}.tag.broad -g hs --nomodel --shift -75 --extsize 150 --keep-dup all --broad --broad-cutoff 0.1
*
*    /home/asd2007/ATACseq/broadpeak.py  ${Sample}.tag.broad_peaks.broadPeak ${Sample}.tn5.broadPeak 
*    """
*        }
*
*
*process signalTrack {
*    // tag "$Sample"
*
*    cpus 16
*
*        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false
*
*        input:
*        set Sample, file(sbam) from bamforsignal
*        val sz from sizefactors
*
*
*        output:
*        set Sample, file("${Sample}.sizefactors.bw") into signal
*
*        script:
*        """
*        samtools index ${sbam} &&
*        bamCoverage --bam ${sbam} --binSize 5 \
*            --outFileFormat bigwig --smoothLength 100 \
*            --normalizeUsingRPKM --scaleFactor ${sz} \
*            --maxFragmentLength 150 \
*            -o ${Sample}.sizefactors.bw --centerReads --extendReads --numberOfProcessors ${task.cpus}
*
*        """
*        }
*/



workflow.onComplete {
    def subject = 'pipeline execution'
    def recipient = 'ashley.doane@gmail.com'

    ['mail', '-s', subject, recipient].execute() << """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}
