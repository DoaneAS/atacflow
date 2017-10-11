<!--# ![atacflow](http://physiology.med.cornell.edu/faculty/elemento/lab/wcmc_logo.gif) -->
# AtacFlow
## Analysis pipeline for ATAC-seq data using Nextflow

This pipeline inspired by and based on the [ENCODE ATAC-seq processubg pipeline](https://www.encodeproject.org/atac-seq/) and
the *prototype* ATAC-seq pipeline
developed by [Anshul Kundaje's lab](https://github.com/kundajelab/atac_dnase_pipelines) at Stanford University

# Installation
* Install [Nextflow](https://www.nextflow.io)
* Clone repository 
  * using nextflow: ```nextflow clone DoaneAS/atacflow ./```
  * or using git: ```git clone https://github.com/DoaneAS/atacflow.git```
* Trimmed read pairs in `./data/Samplename/...trim.fastq.gz`
`python bin/makeIndex.py`

`./callFlow.sh`
