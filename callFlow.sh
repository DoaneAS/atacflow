#!/bin/bash

nextflow run -resume main.nf -with-trace -with-timeline -with-dag flowchart.html
