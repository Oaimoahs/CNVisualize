# Software Requirement Specification

for: CNVisualize

by: Shaomiao Xia

**1. Introduction**

A tool to generate heatmaps of number of reads and the absolute log2 ratio of the two phase states in the aspect of small windows (fixed-length windows, 100-SNP windows, or the phase state related windows) for single-cell RNA-seq data. These heatmaps would be helpful to find some patterns of CNVs in an intuitive way.

**2. Overall description**

Input: A 10X chromium single-cell RNA-seq bam file, a list of barcodes, and the corresponding vcf file.

Output: A set of heatmaps or more informative plots with different choices of features.

**3. Requirements**

The corresponding ginkgo call VCF file and the barcode list file should be input along with the 10X chromium single-cell RNA-seq bam file.
