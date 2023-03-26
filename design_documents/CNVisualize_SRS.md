# Software Requirement Specification

for: CNVisualize

by: Shaomiao Xia

**1. Introduction**

A tool to generate heatmaps of number of reads in the aspect of small fixed-length windows for single-cell RNA-seq data. These heatmaps would be helpful to find some patterns of CNVs in an intuitive way.

**2. Overall description**

Input: single-cell bam files split from A 10X chromium single-cell RNA-seq bam file, the corresponding bam index file, and the list of barcodes.

Output: A set of heatmaps or more informative plots with different choices of features.

**3. Requirements**

The corresponding bam index file "*.bai" should be include along with its original bam file.

Some third party packages are essential for this tool, including samtools, pandas, numpy, seaborn, scipy, csv, logging, and argparse.
