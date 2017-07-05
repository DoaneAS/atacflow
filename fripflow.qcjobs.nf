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

params.index = 'sindex2.lowpass.tsv'


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

atacs = Channel
       .from(index.readLines())
       .map { line ->
       def list = line.split()
              def bed = file(list[0])
              def peaks = file(list[1])
              def sname = list[2]
              def dprefix = file(list[3])
       println bed
       println peaks
       [ sname, bed, peaks, dprefix ]
}



process frip {

  publishDir "$results_path/frip", mode: 'copy', overwrite: true

  input:
  set sname, file(bed), file(peaks), file(dprefix) from atacs

  output:
  set sname, file("${sname}.frip.txt") into frips

  script:
  """
  getFripQC.py --bed ${bed} --peaks ${peaks} --out ${sname}.frip.txt
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
