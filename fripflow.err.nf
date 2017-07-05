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
params.index = 'samps.txt'

////// Print parameters ///////
log.info ''
log.info ''
log.info 'A T A C - Q C ~ FRiP Scores'
log.info '---------------------------------'
log.info ''
log.info "Index File             : ${params.index}"
log.info ''

       index = file(params.index)

results_path = "$PWD/results"

// Clear pipeline.db file

////// Check input parameters //////

if (!params.index) {
  exit 1, "Please specify the input table file"
}



atacs2 = Channel
    .from(index.readLines())
    .map { line ->
           def list = line.split(',')
           def sname = list[0]
           def path = file(list[1])
           def bed = file("$path/$sname/*.tn5.tagAlign.gz")
           def peaks =file("$path/$sname/*tag.broad_peaks.broadPeak")
              // def readsp = "$path/*{R1,R2}.trim.fastq.gz"
              //  def R1 = file(list[2])
              //    def R2 = file(list[3])
              def message = '[INFO] '
              log.info message
           [ sname, bed, peaks ]
}



process frip {


         input:
         set sname, file(bed), file(peaks) from atacs2

         output:
         set sname, file("${sname}.frip.txt") into frips

         script:
         """
         /home/asd2007/ATACseq/getFrip.py --bed ${bed} --peaks ${peaks} --out ${sname}.frip.txt
         """
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
