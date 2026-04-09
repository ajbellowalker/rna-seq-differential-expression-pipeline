#!/bin/bash

# Align reads using Hisat2
hisat2 -x genome_index -U sample.fastq.gz -S output.sam