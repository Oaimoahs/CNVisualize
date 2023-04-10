# CNVisualize
A tool to plot number of reads in small windows for multiple samples as a heatmap, for single-cell RNA-seq data. These heatmaps would be helpful to find some patterns of CNVs.

For example, the heatmap will indicate the number of reads in each small binned windows in the whole genome, then it would be clear to see which windows contain significantly less reads, and then it might indicate there there is a potential CNV in that region.

The tool takes split single-cell RNA-seq bam files, a list of barcodes, reference genome chromosome sizes, and the window size for binnig as input, and the output would be a heatmap contains the number of reads values for each cell sample.

DOWNLOAD the split bam files in this link below:
https://www.dropbox.com/sh/ke5s2w7gb96m09g/AACM1vP24MTBQDFTNwUgjlWJa?dl=0
