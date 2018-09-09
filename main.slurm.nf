#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
         ELEMENTO - MELNICK LABS ||  ATAC-SEQ   BEST   PRACTICE
========================================================================================
#### Authors
Ashley Stephen Doane <asd2007@med.cornell.edu>

*/
//params.index = '/home/asd2007/Scripts/nf/fripflow/sindex.tsv'



version = 1.2

// SET PARAMS

params.name = 'ecadd4547'
params.index = 'sampleIndex.csv'
params.outdir = "$baseDir/results"
params.name = false
params.project = false
       //params.genome = 'hg38'
params.genomes = []
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.black = params.genome ? params.genomes[ params.genome ].BLACK ?: false : false
params.notrim = false
params.saveReference = true
params.saveTrimmed = false
params.saveAlignedIntermediates = false
params.broad = false
params.outdir = './results'
params.email = 'ashley.doane@gmail.com'
params.lncaprefpeak = "$baseDir/data/lncapPeak.narrowPeak" 
params.bcellrefpeak = "$baseDir/data/gcb.tn5.broadPeak" 
params.bt2_index = params.genome ? params.genomes[ params.genome ].bt2 ?: false : false
params.species = params.genome ? params.genomes[ params.genome ].spec ?: false : false
params.picardconfig = params.genome ? params.genomes[ params.genome ].picardconf ?: false : false
params.DNASE = params.genome ? params.genomes[ params.genome ].DNASE ?: false : false
params.TSS_ENRICH = params.genome ? params.genomes[ params.genome ].TSS_ENRICH ?: false : false
params.PROM = params.genome ? params.genomes[ params.genome ].PROM ?: false : false
params.ENH = params.genome ? params.genomes[ params.genome ].ENH ?: false : false
params.REG2MAP = params.genome ? params.genomes[ params.genome ].REG2MAP ?: false : false
params.ROADMAP_META = params.genome ? params.genomes[ params.genome ].ROADMAP_META ?: false : false
params.REF = params.genome ? params.genomes[ params.genome ].REF ?: false : false
params.REG2MAP_BED = params.genome ? params.genomes[ params.genome ].REG2MAP_BED ?: false : false
params.chrsz = params.genome ? params.genomes[ params.genome ].chrsz ?: false : false
params.encodedhs = params.genome ? params.genomes[ params.genome ].ENCODEDHS ?: false : false
       //params.picardconfig="/athena/elementolab/scratch/asd2007/reference/hg38/picardmetrics.conf"


if( params.bwa_index ){
    bwa_index = Channel
        .fromPath(params.bwa_index)
        .ifEmpty { exit 1, "BWA index not found: ${params.bwa_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
        if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
                              } else {
    exit 1, "No reference genome specified!"
}


if( params.bt2_index ){
    bt2_index = Channel
        .fromPath(params.bt2_index)
        .ifEmpty { exit 1, "Bowtie2 index not found: ${params.bt2_index}" }
} else if ( params.fasta ){
    fasta = file(params.fasta)
        if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
                                  } else {
    exit 1, "No reference genome specified!"
        }




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

       //bwaref = file(params.ref)
results_path = file(params.outdir)
picardconf = file(params.picardconfig)
blacklist = file(params.blacklist)
bcellrefpeaks = file(params.bcellrefpeak)
dnase = file(params.DNASE)
tssenrich = file(params.TSS_ENRICH)
prom = file(params.PROM)
enh = file(params.ENH)
reg2map = file(params.REG2MAP)
reg2mapbed = file(params.REG2MAP_BED)
roadmapmeta = file(params.ROADMAP_META)
ref = file(params.REF)
blackqc = file(params.black)
black = file(params.black)
encodedhs = file(params.encodedhs)
chrsz = file(params.chrsz)
species = params.species


 /*
 * PREPROCESSING - Build BWA index
 */
if(!params.bwa_index && fasta){
    process makeBWAindex {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        """
        spack load bwa
        mkdir BWAIndex
        bwa index -a bwtsw $fasta
        """
    }
}



/*
 * PREPROCESSING - Build Bt2 index
 */
if(!params.bt2_index && fasta){
    process makeBT2index {
        tag fasta
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
        saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta

        output:
        file "BWAIndex" into bwa_index

        script:
        """
        spack load bowtie2
        mkdir BT2Index
        bowtie2-build $fasta genome
        """
    }
}

////// Check input parameters //////

if (!params.index) {
  exit 1, "Please specify the input table file"
}


sizefactors = Channel.from(1)





fastq = Channel
.from(index.readLines())
.map { line ->
       def list = line.split(',')
       def Sample = list[0]
       def path = file(list[1])
       //def reads = file("$path/*_{R1,R2}_001.fastq.gz")
       def reads = file("$path/*.{R1,R2}.fastq.gz")
       // def readsp = "$path/*{R1,R2}.trim.fastq.gz"
       //  def R1 = file(list[2])
       //    def R2 = file(list[3])
       def message = '[INFO] '
       log.info message
       [ Sample, path, reads ]
}



/*
 * STEP 2 - Trim Galore!
 */
//if(params.notrim){
///    trimmed_reads = read_files_trimming
//} else {
process ngtrim {
    tag "$Sample"
    publishDir "$results_path/$Sample/$Sample", mode: 'copy'

    cpus 4
    executor 'slurm'
        ////clusterOptions '-l h_vmem=2G -l h_rt=24:00:00 -l athena=true'
        //scratch true
    time '4h'
    memory '2 GB'

    input:
    set Sample, file(path), file(reads) from fastq

    output:
        set Sample, file(path), file('*_{1,2}.fastq.gz') into trimmed_reads
        //file '*trimming_report.txt' into trimgalore_results
        //file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

    script:
        def R1 = reads[0]
        def R2 = reads[1]
        """
        NGmerge -z -n ${task.cpus} -a -1 $R1 -2 $R2 -o ${Sample}
        """
            }
//}






process bt2 {
    tag "$Sample"
        publishDir "$results_path/$Sample/$Sample", mode: 'copy'
        //conda 'bowtie2 samtools'

        cpus 8
        executor 'slurm'
        memory '4 GB'
        //penv 'smp'
        //clusterOptions '-l h_vmem=4G -l h_rt=24:00:00 -l athena=true -R y'
        scratch true

        input:
        set Sample, file(path), file(reads) from trimmed_reads
        file index from bt2_index.first()
        //file bwaref from bwa_index.collect()
        //file(bwaref) from bwaref

        output:
        set Sample, file("${Sample}.bam") into newbam
        set Sample, file("${Sample}.bt2.log") into bt2log
   

        script:
        def R1 = reads[0]
        def R2 = reads[1]
        """
        #!/bin/bash -l
        set -o pipefail
        spack load bowtie2@2.3.4.1
        spack load samtools
        bowtie2 -X2000 -x ${index}/genome --local -p 7 -1 ${R1} -2 ${R2} 2> ${Sample}.bt2.log | samtools view -bS -q 30 - > ${Sample}.bam
        """
        }



process processbam {
    tag "$Sample"
    publishDir "$results_path/$Sample/$Sample", mode: 'copy'

    executor 'slurm'
    cpus 8
    time '8h'
    memory '4 GB'
    //penv 'smp'
    //clusterOptions '-l h_vmem=4G -l h_rt=16:00:00 -l athena=true'
    scratch true
    // cpus 8


    input:
    set Sample, file(nbam) from newbam
    file(BLACK) from blacklist

    output:
    set Sample, file("${Sample}.sorted.bam") into sortedbam
    set Sample, file("${Sample}.sorted.bam") into sortbamqc
    set Sample, file("${Sample}.sorted.bam") into sortedbamqc
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbam
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into finalbamforqc
    set Sample, file("${Sample}.sorted.nodup.noM.black.bam") into bamforsignal
    set Sample, file("${Sample}*.pbc.qc") into pbcqc
    set Sample, file("${Sample}*.dup.qc") into dupqc
    file("*nsort.fixmate.bam") into fixmatebam
        // file("*window500.hist_data") into hist_data
    file("*window500.hist_graph.pdf") into fragsizes
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbam
        // set Sample, file("${Sample}.sorted.nodup.noM.black.bam"), file("${Sample}.sorted.nodup.noM.black.bam.bai") into bamforsignal
    set Sample, file("${Sample}.nsorted.nodup.noM.bam") into nsortedbamforqc

    script:
    """
    #!/bin/bash -l
    set -o pipefail
    ##. ~/.spackloads.sh
    source activate atacFlow
    processAlignment.nf.sh ${nbam} ${BLACK} 8
    """
}


//hist_data.subscribe { println "Received: " + file(hist_data)}

//fragsizes.subscribe { println "Received: " + file(fragsizes)}



process bam2bed {
    tag "$Sample"
    publishDir  "$results_path/$Sample/$Sample", mode: 'copy'

    input:
    set Sample, file(nsbam) from nsortedbam

    output:
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbedqc
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz") into finalbed
    set Sample, file("${Sample}.nodup.bedpe.gz") into finalbedpe
    set Sample, file("${Sample}.nodup.tn5.tagAlign.gz"), file("${Sample}.nodup.bedpe.gz") into finalbedmacs

    script:
    """
    samtools fixmate ${nsbam} ${Sample}.nsorted.fixmate.nodup.noM.bam
    convertBAMtoBED.sh ${Sample}.nsorted.fixmate.nodup.noM.bam
    cp ${Sample}.nsorted.fixmate.nodup.noM.tn5.tagAlign.gz  ${Sample}.nodup.tn5.tagAlign.gz
    cp ${Sample}.nsorted.fixmate.nodup.noM.bedpe.gz ${Sample}.nodup.bedpe.gz
    """
        }




process callpeaks {
    tag "$Sample"
    publishDir  "$results_path/$Sample/$Sample", mode: 'copy'

    input:
    set Sample, file(rbed), file(rbedpe) from finalbedmacs
    val sp from species

    output:

    set Sample, file("${Sample}*.narrowPeak") into narrowpeakfile
    set Sample, file("${Sample}.tn5.narrowPeak.gz") into narrowpeak
    set Sample, file("${Sample}.tag.narrow_summits.bed") into summits
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeak
    set Sample, file("${Sample}.tn5.broadPeak.gz") into broadpeakqc
        // set Sample, file("${Sample}/pseudoreps/${Sample}.tn5.pooled.pf.pval0.1.500K.naive_overlap.narrowPeak.gz") into pooledpeaks

    script:
    """
    #!/bin/bash -l
    source activate atacFlow
    callbedpepeaks.sh ${rbed} ${Sample} ${sp}

    """
        }


process signalTrack {
    tag "$Sample"

    publishDir "$results_path/$Sample/$Sample", mode: 'copy'

        //conda 'bioconda::deeptools=3.1.2'
    executor 'slurm'
    //penv 'smp'
    memory { 3.GB * task.attempt }
    //clusterOptions '-l h_rt=48:00:00 -l athena=true'
    scratch true
    cpus 18
    errorStrategy { task.exitStatus == 140 ? 'retry' : 'terminate' }

    input:
    set Sample, file(sbam) from bamforsignal
        //val sz from sizefactors


    output:
    set Sample, file("${Sample}.bpm.sizefactors.bw") into insertionTrackbw

    script:
    """
    #!/bin/bash -l
    source activate deep

    getbamcov.sh ${sbam} ${Sample} ${task.cpus}

    source deactivate
    """
    }


finalbedpe.mix(broadpeak)
    .groupTuple(sort: true)
    .set{ fripin }

process frip {
    tag "$Sample"

    publishDir "$results_path/$Sample/frip", mode: 'copy'

    input:
    set Sample, file(file_list) from fripin
        //set Sample, file(peaks) from broadpeak
        //   file(lncapref) from lncaprefpeaks
    file(bcellref) from bcellrefpeaks
    file(encodedhs) from encodedhs

    output:
    set Sample, file("${Sample}.frip.txt") into frips
    set Sample, file("${Sample}.bcellref.frip.txt") into frips2
    set Sample, file("${Sample}.encodedhs.frip.txt") into frips3

    script:
    """
    #!/bin/bash -l

    source activate atacFlow

    getFripQC.py \\
    --bed ${Sample}.nodup.bedpe.gz --peaks ${Sample}.tn5.broadPeak.gz --out ${Sample}.frip.txt


    getFripQC.py \\
    --bed ${Sample}.nodup.bedpe.gz --peaks ${bcellref} --out ${Sample}.bcellref.frip.txt


    getFripQC.py \\
    --bed ${Sample}.nodup.bedpe.gz --peaks ${encodedhs} --out ${Sample}.encodedhs.frip.txt

    source deactivate

        """
}





process picardqc {
    tag "$Sample"
    publishDir "$results_path/$Sample/qc", mode: 'copy'

    input:
    set Sample, file(sortbamqc) from sortedbamqc
    file(picardconfig) from picardconf 

    output:
    set Sample, file("QCmetrics/${Sample}.picardcomplexity.qc") into picardcomplexity
        //set Sample, file(sortbamqc) into sortbamqc

    script:
    """
    mkdir -p QCmetrics
    source /home/asd2007/Scripts/picard.env 
    spack load jdk@8u172-b11
    spack load samtools
    picard EstimateLibraryComplexity I=${sortbamqc} O=${Sample}.EstimateLibraryComplexity.log
    cp *.EstimateLibraryComplexity.log QCmetrics/${Sample}.picardcomplexity.qc
    """

}


finalbamforqc.mix(nsortedbamforqc)
    .mix(broadpeakqc)
    .mix(finalbedqc)
    .mix(sortbamqc)
    .mix(insertionTrackbw)
    .mix(picardcomplexity)
    .mix(pbcqc)
    .mix(dupqc)
    .mix(frips)
    .mix(bt2log)
    .groupTuple(sort: true)
    .view()
    .set{ qcin }


process atacqc {
    tag "$Sample"

    publishDir "$results_path/$Sample/qc", mode: 'copy'
        //conda 'bds_atac_requirements.txt'

    cpus 4

    input:
    set Sample, file(file_list) from qcin
    file(dnase) from dnase
    file(tssenrich) from tssenrich
    file(prom) from prom
    file(enh) from enh
    file(reg2map) from reg2map
    file(roadmapmeta) from roadmapmeta
    file(ref) from ref
    file(black) from blackqc
    file(reg2mapbed) from reg2mapbed

    output:
    set Sample, file("${Sample}*.preseq.log"), file("${Sample}*_qc.txt"), file("${Sample}*large_vplot.png"), file("${Sample}*vplot.png"), file("${Sample}_qc.trad.txt"), file("${Sample}*qc.html"), file("*qc.save") into qcdat
    set Sample, file("*.log"), file("*qc") into logs


    script:
    """
    #!/bin/bash -l

    OUTPREFIX=${Sample}
    INPREFIX=${Sample}
    SAMPLE=${Sample}
    PBC=${Sample}.pbc.qc
    echo ${Sample} ${file_list}
    ##spack load jdk

    spack load samtools
    samtools index -@4 ${Sample}.sorted.nodup.noM.black.bam
    samtools index -@4 ${Sample}.sorted.bam


    source activate atacFlow && ${baseDir}/bin/run_ataqc.athena.py --workdir \$PWD  \\
    --outdir \$PWD \\
    --outprefix ${Sample} \\
    --genome hg38 \\
    --ref ${ref} --tss ${tssenrich} \\
    --dnase ${dnase} \\
    --blacklist ${black} \\
    --prom ${prom} \\
    --enh ${enh} \\
    --reg2map ${reg2map} \\
    --reg2map_bed ${reg2mapbed} \\
    --meta ${roadmapmeta} \\
    --alignedbam ${Sample}.sorted.bam  \\
    --alignmentlog ${Sample}.align.log \\
    --coordsortbam ${Sample}.sorted.bam \\
    --duplog ${Sample}.dup.qc \\
    --pbc ${Sample}.pbc.qc \\
    --finalbam ${Sample}.sorted.nodup.noM.black.bam \\
    --finalbed ${Sample}.nodup.tn5.tagAlign.gz \\
    --bigwig ${Sample}.bpm.sizefactors.bw \\
    --peaks ${Sample}.tn5.broadPeak.gz \\
    --naive_overlap_peaks ${Sample}.tn5.broadPeak.gz \\
    --idr_peaks ${Sample}.tn5.broadPeak.gz  --processes ${task.cpus}

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

