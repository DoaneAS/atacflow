#!/bin/bash

nextflow run -with-trace -with-timeline -with-dag flowchart.html \
         main.nf --index sampleIndex.csv --genome hg38
