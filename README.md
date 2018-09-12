<!--# ![atacflow](http://physiology.med.cornell.edu/faculty/elemento/lab/wcmc_logo.gif) -->
# AtacFlow
## Analysis pipeline for ATAC-seq data using Nextflow

This pipeline inspired by and based on the [ENCODE ATAC-seq processubg pipeline](https://www.encodeproject.org/atac-seq/) and
the *prototype* ATAC-seq pipeline
developed by [Anshul Kundaje's lab](https://github.com/kundajelab/atac_dnase_pipelines) at Stanford University

## Installation
* Install [Nextflow](https://www.nextflow.io)
* Clone repository 
  * using nextflow: ```nextflow clone DoaneAS/atacflow ./```
  * or using git: ```git clone https://github.com/DoaneAS/atacflow.git```
* Install conda dependencies:
   ```
   conda update conda
   conda env create --file requirements.atacFlow.yml
   conda env create --file deep.yml
   ```

## Setup data
* ATAC-seq reads go in ```data/<Sample>/*_001.fastq.gz```
  * Concatenate read pairs per sample ```parallel -j8 './bin/catlanes.sh {}' ::: data/Sample*```
* Create sample index: `python bin/makeIndex.py`

## Execution  
```
nextflow run -with-trace -with-dag flow.html main.nf --index sampleIndex.csv --genome hg38
```  
* supported genomes on panda WCM cluster:  hg38, mm10
