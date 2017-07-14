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

params.name = 'ecadlowpass'
       //params.index = 'sampleIndex.csv'
   //params.index = 'indexpooled.tsv'
       //params.index = 'Sample_Ly7_pooled_100k'
params.index = "$baseDir/indexTest.csv"
       //params.index = 'sampleIndexjc.csv'
params.genome = 'hg38'
params.blacklist = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
       //genome = file(params.genome)
index = file(params.index)
params.chrsz = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"
params.ref = "/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.project = 'ecadlowpass'
       //params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
       //params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.notrim = false
params.saveReference = false
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.broad = false
params.outdir = './results'
params.email = 'ashley.doane@gmail.com'
params.chromsizes = "/athena/elementolab/scratch/asd2007/reference/hg38/hg38.chrom.sizes"
params.lncaprefpeak = "$baseDir/data/lncapPeak.narrowPeak" 
params.bcellrefpeak = "$baseDir/data/gcb.tn5.broadPeak" 

params.DNASE_BED="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
params.BLACK="/athena/elementolab/scratch/asd2007/reference/hg38/hg38.blacklist.bed.gz"
params.PICARDCONF="/athena/elementolab/scratch/asd2007/reference/hg38/picardmetrics.conf"
params.REF="/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.REFbt2="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.RG="hg38"
params.SPEC="hs"
params.REFGen="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.seq="/athena/elementolab/scratch/asd2007/reference/hg38/seq"
params.gensz="hs"
params.bwt2_idx="/athena/elementolab/scratch/asd2007/reference/hg38/bowtie2_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.REF_FASTA="/athena/elementolab/scratch/asd2007/reference/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
params.species_browser='hg38'
params.TSS_ENRICH="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
params.DNASE='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz'
params.PROM='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz'
params.ENH='athena/elementolab/scratch/asd2007/reference/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz'
params.REG2MAP='/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz'
params.REG2MAP_BED="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
params.ROADMAP_META="/athena/elementolab/scratch/asd2007/reference/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"





custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// Header log info
log.info "=================================================================="
log.info " Elemento-Melnick-Labs-ATACseq: ATAC-Seq Best Practice v${version}"
log.info "==================================================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
       //summary['Reads']        = params.reads
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

bcellrefpeaks = Channel.fromPath(params.bcellrefpeak)

lncaprefpeaks = Channel.fromPath(params.lncaprefpeak)

//    Channel
//       .fromFilePairs( params.reads )                                             
//       .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
//        .set { read_pairs }





process bwamem {

    executor 'sge'
    scratch true
    clusterOptions '-l h_vmem=4G -pe smp 8 -l h_rt=16:00:00 -l athena=true'
    publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false


    input:
    set Sample, file(path), file(reads) from fastq
        //file(bwaref) from bwaref

    output:
    set Sample, file("${Sample}.bam") into newbam

    script:
    """
    #!/bin/bash -l
    set -o pipefail
    spack load bwa
    bwa mem -t 8 -M /athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta $reads | samtools view -bS - > ${Sample}.bam
    """

    }






process processbam {


    // executor 'sge'
        //    clusterOptions '-l h_vmem=4G -pe smp 6 -l h_rt=16:00:00 -l athena=true'
    // scratch true
    cpus 8

    publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false

    input:
    set Sample, file(nbam) from newbam
    file(BLACK) from blacklist

    output:
    set Sample, file("${Sample}.sorted.bam") into sortedbam
    set Sample, file("${Sample}.sorted.bam") into sortedbamqc
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbam
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbamforqc
    file("*.pbc.qc") into pbcqc
    file("*.dup.qc") into dupqc
    file("*nsort.fixmate.bam") into fixmatebam
    file("*window500.hist_data") into hist_data
    file("*window500.hist_graph.pdf") into fragsizes
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbamforqc

    script:
    """
    processAlignment.nf.sh ${nbam} ${BLACK} 8
    ##sambamba sort --memory-limit 38GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbam}
    ##samtools index ${Sample}.sorted.nodup.noM.black.bam
    """
}


hist_data.subscribe { println "Received: " + file(hist_data)}

fragsizes.subscribe { println "Received: " + file(fragsizes)}

/*
*process nsortbam {
*    // tag "$Sample"
*
*        cpus 4
*        memory '20 GB'
*
*        publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false
*
*        input:
*        set Sample, file(finalbam) from finalbam
*
*        output:
*        set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
*        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
*        set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
*        set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbamforqc
*            //val sf into sizefactors
*
*        script:
*            // def sf = 1
*            //def fbam = file(finalbam)
*            // finalbam.renameTo("${Sample}.sorted.nodup.noM.black.bam")
*
*        """
*        samtools index ${finalbam}
*        sambamba sort --memory-limit 18GB -n -t ${task.cpus} --out ${Sample}.nsorted.nodup.noM.bam ${finalbam}
*        samtools index ${Sample}.sorted.nodup.noM.black.bam
*        """
*}
*/




process bam2bed {

    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false

    input:
    set Sample, file(nsbam) from nsortedbam

    output:
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbedqc
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

    publishDir  "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false

    input:
    set Sample, file(rbed) from finalbed

    output:
    set Sample, file("${Sample}.tn5.narrowPeak.gz") into narrowpeak
    set Sample, file("${Sample}.tag.narrow_summits.bed") into summits
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeak
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeakqc

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

    publishDir "$results_path/$Sample/$Sample", mode: 'copy', overwrite: false
        // executor 'sge'
        //clusterOptions '-l h_vmem=5G -pe smp 8 -l h_rt=16:00:00 -l athena=true'
        //scratch true
    cpus 8

    input:
    set Sample, file(sbam) from bamforsignal
        //val sz from sizefactors


    output:
    set Sample, file("*.bw") into insertionTrack

    script:
    """
    spack load samtools
    samtools index ${sbam}
    bamCoverage --bam ${sbam} --binSize 20 --outFileFormat bigwig --smoothLength 120 \
        --normalizeUsingRPKM \
        --maxFragmentLength 150 \
        -o ${Sample}.sizefactors.bw --centerReads --extendReads --numberOfProcessors 8

    """
    }




process frip {

    publishDir "$results_path/$Sample/qc", mode: 'copy', overwrite: true

        input:
        set sname, file(bed) from finalbedpe
        set Sample, file(peaks) from broadpeak
        file(lncapref) from lncaprefpeaks
        file(bcellref) from bcellrefpeaks

        output:
        set sname, file("${sname}.frip.txt") into frips
        set sname, file("${sname}.lncapref.frip.txt") into frips2
        set sname, file("${sname}.bcellref.frip.txt") into frips3

        script:
        """
        getFripQC.py --bed ${bed} --peaks ${peaks} --out ${sname}.frip.txt

        getFripQC.py --bed ${bed} --peaks ${lncapref} --out ${sname}.lncapref.frip.txt

        getFripQC.py --bed ${bed} --peaks ${bcellref} --out ${sname}.bcellref.frip.txt
        """
        }





process picardqc {

    publishDir "$results_path/$Sample/qc", mode: 'copy', overwrite: true

    input:
    set Sample, file(sortbamqc) from sortedbamqc


    output:
    set Sample, file("QCmetrics/${Sample}.picardcomplexity.qc") into picardcomplexity
    set Sample, file(sortbamqc) into sortbamqc

    script:
    """
    mkdir -p QCmetrics
    picardmetrics run -f $PICARDCONF -o QCmetrics $sortbamqc
    cp QCmetrics/*.EstimateLibraryComplexity.log QCmetrics/${Sample}.picardcomplexity.qc
    """

}



 
process atacqc {
 
     publishDir "$results_path/$Sample/qc", mode: 'copy', overwrite: true
 
         input:
         file(pbc) from pbcqc
         set Sample, file(finalbamqc) from finalbamforqc
         set Sample, file(nbamforqc) from nsortedbamforqc
         set Sample, file(broadpeaks) from broadpeakqc
         set Sample, file(finalbedqc) from finalbedqc
         set Sample, file(sortbamqc) from sortbamqc
         set Sample, file(insertioTrack) from insertionTrack
         file(dupqc) from dupqc


         output:
         set Sample, file("$Sample*.preseq.dat"), file("$Sample*_qc.txt"), file("$Sample*large_vplot.png"), file("$Sample*vplot.png") into qcdat1
         set Sample, file("*.log"), file("*qc") into logs


         script:
         """
         samtools index $finalbamqc
         samtools index $sortbamqc
         python run_ataqc.athena.py --workdir \${PWD} \
             --outdir qc \
             --outprefix ${Sample} \
             --genome \${GENOME} \
             --ref \${REF} --tss \$TSS_ENRICH \
             --dnase \${DNASE} \
             --blacklist ${BLACK} \
             --prom \$PROM \
             --enh \${ENH} \
            --reg2map \${REG2MAP} \
            --meta \${ROADMAP_META} \
            --alignedbam $sortbamqc  \
            --coordsortbam $sortbamqc \
            --duplog $dupqc \
            --pbc $pbc \
            --finalbam $finalbamqc \
            --finalbed $finalbedqc \
            --bigwig $insertionTrack \
            --peaks $broadpeaks \
            --naive_overlap_peaks $broadpeaks \
            --idr_peaks $broadpeaks --processes 4

         """

         }



workflow.onComplete {
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
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
 *
 *workflow.onComplete {
 *
 *    // Set up the e-mail variables
 *    def subject = "[OElab-ATACseq] Successful: $workflow.runName"
 *    if(!workflow.success){
 *      subject = "[OElab-ATACseq] FAILED: $workflow.runName"
 *    }
 *    def email_fields = [:]
 *    email_fields['version'] = version
 *    email_fields['runName'] = custom_runName ?: workflow.runName
 *    email_fields['success'] = workflow.success
 *    email_fields['dateComplete'] = workflow.complete
 *    email_fields['duration'] = workflow.duration
 *    email_fields['exitStatus'] = workflow.exitStatus
 *    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
 *    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
 *    email_fields['commandLine'] = workflow.commandLine
 *    email_fields['projectDir'] = workflow.projectDir
 *    email_fields['summary'] = summary
 *    email_fields['summary']['Date Started'] = workflow.start
 *    email_fields['summary']['Date Completed'] = workflow.complete
 *    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
 *    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
 *    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp
 *    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
 *    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
 *    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
 *    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
 *    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
 *    if(workflow.container) email_fields['summary']['Singularity image'] = workflow.container
 *
 *    // Render the TXT template
 *    def engine = new groovy.text.GStringTemplateEngine()
 *    def tf = new File("$baseDir/assets/email_template.txt")
 *    def txt_template = engine.createTemplate(tf).make(email_fields)
 *    def email_txt = txt_template.toString()
 *
 *    // Render the HTML template
 *    def hf = new File("$baseDir/assets/email_template.html")
 *    def html_template = engine.createTemplate(hf).make(email_fields)
 *    def email_html = html_template.toString()
 *
 *    // Render the sendmail template
 *    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
 *    def sf = new File("$baseDir/assets/sendmail_template.txt")
 *    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
 *    def sendmail_html = sendmail_template.toString()
 *
 *    // Send the HTML e-mail
 *    if (params.email) {
 *        try {
 *          // Try to send HTML e-mail using sendmail
 *          [ 'sendmail', '-t' ].execute() << sendmail_html
 *          log.debug "[NGI-ChIPseq] Sent summary e-mail using sendmail"
 *        } catch (all) {
 *          // Catch failures and try with plaintext
 *          [ 'mail', '-s', subject, params.email ].execute() << email_txt
 *          log.debug "[OElab-ATACseq] Sendmail failed, failing back to sending summary e-mail using mail"
 *        }
 *        log.info "[OElab-ATACseq] Sent summary e-mail to $params.email"
 *    }
 *
 *    // Write summary e-mail HTML to a file
 *    def output_d = new File( "${params.outdir}/Documentation/" )
 *    if( !output_d.exists() ) {
 *      output_d.mkdirs()
 *    }
 *    def output_hf = new File( output_d, "pipeline_report.html" )
 *    output_hf.withWriter { w -> w << email_html }
 *    def output_tf = new File( output_d, "pipeline_report.txt" )
 *    output_tf.withWriter { w -> w << email_txt }
 *
 *    log.info "[OElab-ATACseq] Pipeline Complete"
 *}
 */

