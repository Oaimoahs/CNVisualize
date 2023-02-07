# CNVisualize
A tool to generate heatmaps of number of reads and the absolute log2 ratio of the two phase states in the aspect of small windows (fixed-length windows, 100-SNP windows, or the phase state related windows) for single-cell RNA-seq data. These heatmaps would be helpful to find some patterns of CNVs.

The tool takes a 10X chromium single-cell RNA-seq bam file, a list of barcodes, and the corresponding vcf file as input, and the output would be a set of heatmaps or more informative plots.
