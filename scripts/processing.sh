#!/bin/bash

# Convert SAM to BAM
samtools view -bS output.sam > output.bam

# Sort BAM
samtools sort output.bam -o output_sorted.bam

# Index BAM
samtools index output_sorted.bam