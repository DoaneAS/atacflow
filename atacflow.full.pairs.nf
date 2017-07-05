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



version = 1.1

// SET PARAMS 
params.index = 'sampleIndex.csv'
   //params.index = 'indexpooled.tsv'
params.genome = 'hg38'
params.blacklist = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
genome = file(params.genome)
index = file(params.index)
params.chrsz = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"
params.ref = "/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.name = false
params.project = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.notrim = false
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.broad = false
params.outdir = './results'
params.email = 'ashley.doane@gmail.com'
params.chromsizes = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"


// Header log info
log.info "=================================================================="
log.info " Elemento-Melnick-Labs-ATACseq: ATAC-Seq Best Practice v${version}"
log.info "==================================================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Genome']       = params.genome
if(params.bwa_index)  summary['BWA Index'] = params.bwa_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Save Reference'] = params.saveReference
summary['Save Trimmed']   = params.saveTrimmed
summary['Save Intermeds'] = params.saveAlignedIntermediates
if(params.notrim)       summary['Trimming Step'] = 'Skipped'
if( params.clip_r1 > 0) summary['Trim R1'] = params.clip_r1
if( params.clip_r2 > 0) summary['Trim R2'] = params.clip_r2
if( params.three_prime_clip_r1 > 0) summary["Trim 3' R1"] = params.three_prime_clip_r1
if( params.three_prime_clip_r2 > 0) summary["Trim 3' R2"] = params.three_prime_clip_r2
if(params.email) summary['E-mail Address'] = params.email
if(workflow.commitId) summary['Pipeline Commit']= workflow.commitId
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "===================================="

index = file(params.index)

bwaref = file(params.ref)
results_path = "./results"

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
              def list = line.split(',')
              def Sample = list[0]
              def path = file(list[1])
              def reads = file("$path/*{R1,R2}.trim.fastq.gz")
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
        scratch true
        clusterOptions '-l h_vmem=5G -pe smp 4 -l h_rt=96:00:00 -l athena=true'
            //     publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true


        input:
        set Sample, file(path), file(reads) from fastq
        file(bwaref) from bwaref

        output:
        set Sample, file("${Sample}.bam") into newbam

        script:
        """
        #!/bin/bash -l
        set -o pipefail
        spack load bwa
        bwa mem -t 4 -M /athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta $reads | samtools view -bS - > ${Sample}.bam
        """

    }






process processbam {


    executor 'sge'
    clusterOptions '-l h_vmem=5G -pe smp 4 -l h_rt=16:00:00 -l athena=true'
    scratch true

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

    script:
    """
    #!/bin/bash -l
    processAlignment.nf.sh $nbam $BLACK 4
    """
}



hist_data.subscribe { println "Received: " + file(hist_data)}

fragsizes.subscribe { println "Received: " + file(fragsizes)}


process nsortbam {
    // tag "$Sample"

        cpus 4
        memory '20 GB'

        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true

        input:
        set Sample, file(finalbam) from finalbam

        output:
        set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
        set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
        set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbamforqc
            //val sf into sizefactors

        script:
            // def sf = 1
            //def fbam = file(finalbam)
            // finalbam.renameTo("${Sample}.sorted.nodup.noM.black.bam")

        """
        samtools index ${finalbam}
        sambamba sort --memory-limit 18GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbam}

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

        executor 'sge'
        clusterOptions '-l h_vmem=5G -pe smp 4 -l h_rt=16:00:00 -l athena=true'
        scratch true
        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: true


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
            -o ${Sample}.sizefactors.bw --centerReads --extendReads --numberOfProcessors 4

        """
        }




process frip {

    publishDir "$results_path/$Sample/qc", mode: 'copy', overwrite: true

        input:
        set sname, file(bed) from finalbedpe
        set Sample, file(peaks) from broadpeak

        output:
        set sname, file("${sname}.frip.txt") into frips

        script:
        """
        getFripQC.py --bed ${bed} --peaks ${peaks} --out ${sname}.frip.txt
        """
        }



process atacqc {

    publishDir "$results_path/$Sample/qc", mode: 'copy', overwrite: true

        input:
        file(pbcqc)
        set Sample, file(finalbamforqc) from nsortbam


        output:
        set Sample, file("QCmetrics/${Sample}.picardcomplexity.qc") into picardcomplexity.qc
        set Sample, file("*.preseq.dat"), file("*_qc.txt") into atacqc

        script:
        """
        ./run_atacqc.athena.nf.sh -s ${Sample} -g hg38
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




/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[OElab-ATACseq] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[OElab-ATACseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Singularity image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.debug "[NGI-ChIPseq] Sent summary e-mail using sendmail"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.debug "[OElab-ATACseq] Sendmail failed, failing back to sending summary e-mail using mail"
        }
        log.info "[OElab-ATACseq] Sent summary e-mail to $params.email"
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[OElab-ATACseq] Pipeline Complete"
}

