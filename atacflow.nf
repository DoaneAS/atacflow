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

params.dbFile = 'chipseq-pipeline.db'
params.genome = '/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
params.genomeIndex = '/athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta'
params.genomeSize = 'hs'
       //params.fragmentLength = 200
params.help = false
params.index = 'samples2.tsv'
params.minMatchedBases = 0.8
params.mismatches = 2
params.multimaps = 10
params.qualityThreshold = 26
params.rescale = false
params.removeDuplicates = false
params.shift = false

//print usage
if (params.help) {
    log.info ''
    log.info 'C H I P - N F ~ ChIP-seq Pipeline'
    log.info '---------------------------------'
    log.info 'Run ChIP-seq analyses on a set of data.'
    log.info ''
    log.info 'Usage: '
    log.info '    chipseq-pipeline.nf --index TSV_FILE --genome GENOME_FILE [OPTION]...'
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --index TSV_FILE                    Tab separted file containing information about the data.'
    log.info '    --genome GENOME_FILE                Reference genome file.'
    log.info '    --genome-index GENOME_INDEX_ FILE   Reference genome index file.'
    log.info '    --genome-size GENOME_SIZE           Reference genome size for MACS2 callpeaks. Must be one of' 
    log.info '                                        MACS2 precomputed sizes: hs, mm, dm, ce. (Default: hs)'
    log.info '    --mismatches MISMATCHES             Sets the maximum number/percentage of mismatches allowed for a read (Default: 2).'
    log.info '    --multimaps MULTIMAPS               Sets the maximum number of mappings allowed for a read (Default: 10).'
    log.info '    --min-matched-bases BASES           Sets the minimum number/percentage of bases that have to match with the reference (Default: 0.80).'
    log.info '    --quality-threshold THRESHOLD       Sets the sequence quality threshold for a base to be considered as low-quality (Default: 26).'
    log.info '    --fragment-length LENGTH            Sets the fragment length globally for all samples (Default: 200).'
    log.info '    --remove-duplicates                 Remove duplicate alignments instead of just flagging them (Default: false).'
    log.info '    --rescale                           Rescale peak scores to conform to the format supported by the'
    log.info '                                        UCSC genome browser (score must be <1000) (Default: false).'
    log.info '    --shift                             Move fragments ends and apply global extsize in peak calling. (Default: false).'
    log.info ''
    exit 1
}
params.dbFile = 'pipeline.db'

genome_file = file(params.genome)

index = file(params.index)

// Clear pipeline.db file
pdb = file(params.dbFile)
pdb.write('')

////// Check input parameters //////
if (!params.genome) {
  exit 1, "Please specify a genome file"
}

if (!params.index) {
  exit 1, "Please specify the input table file"
}
////// End of input parameters check ////////

////// Print parameters ///////
log.info ''
log.info 'C H I P - N F ~ ChIP-seq Pipeline'
log.info '---------------------------------'
log.info ''
log.info "General parameters"
log.info '------------------'
log.info ''
log.info "Index File             : ${params.index}"
log.info "Genome File            : ${params.genome}"
log.info "Genome Index File      : ${params.genomeIndex ?: '-'}"
log.info "MACS2 Genome Size      : ${params.genomeSize}"
log.info "Global Fragment Length : ${params.fragmentLength}"
log.info "Database file          : ${params.dbFile}"
log.info "Remove Duplicates      : ${params.removeDuplicates}"
log.info "Shift                  : ${params.shift}"
log.info "Rescale Peaks          : ${params.rescale}"
log.info ''
log.info "Mapping parameters"
log.info '------------------'
log.info ''
log.info "Max Mismatches         : ${params.mismatches}"
log.info "Max Multimaps          : ${params.multimaps}"
log.info "Minimum Matched Bases  : ${params.minMatchedBases}"
log.info "Low Quality Threshold  : ${params.qualityThreshold}"
log.info ''

genome = file(params.genome)

fastqs = Channel
.from(index.readLines())
.map { line ->
  def list = line.split()
  def mergeId = list[0]
  def id = list[1]
  def read1 = file(list[2])
  def read2 = file(list[3])
       //def controlId = list[3]
       //def mark = list[4]
       //def fragLen = list.size() == 6 ? list[5] as Integer : -1
  def message = '[INFO] '
  [ mergeId, id, read1, read2 ]
}

GenomeIdx = Channel.fromPath(params.genome)


process mapping {

  cpus 2

  input:
  file genomeIndex from GenomeIdx.val
  set mergeId, prefix, file(fastq1), file(fastq2) from fastqs

  output:
  set mergeId, prefix, file("${prefix}.bam") into bams

  script:
  def cpus = task.cpus
  def memory = task.memory
  def readGroup = "ID=${prefix},SM=${mergeId}"
  """
  bwa mem -t ${cpus} -M /athena/elementolab/scratch/asd2007/reference/hg38/bwa_index/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta ${fastq1} ${fastq2} \
    | samtools view -@ ${cpus} -bS - >  ${prefix}.bam
  """
}


process processBam {
    echo true

    input:
    set mergID, prefix file("${prefix}.bam") from bams

    output:
    set mergID, prefix, file("${prefix}.sorted.bam"), file("${prefix}.sorted.nodup.noM.black.bam")i 


}



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
