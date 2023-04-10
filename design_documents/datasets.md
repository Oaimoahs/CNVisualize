# Brief description of the datasets I will use for this project

The dataset I will use for this project is from Ryan Mills' Lab:

1. **Real dataset for answering a biological question using the tool.** The single-cell bams split from a bulk 10X chromium single-cell RNA-seq bam file (20 cells were random selected, data is available here: https://www.dropbox.com/sh/ke5s2w7gb96m09g/AACM1vP24MTBQDFTNwUgjlWJa?dl=0).

2. The corrsponding bam index file (data is available here: https://www.dropbox.com/sh/ke5s2w7gb96m09g/AACM1vP24MTBQDFTNwUgjlWJa?dl=0) .

3. The corresponding barcode list file (For the sample data `/CNVisualize/data/sample_barcodes_test.txt`, for the real dataset for answering a biological question using this tool `CNVisualize/data/sample_barcodes.txt`).

4. The chromosome sizes file of the used reference genome (hg19) (`CNVisualize/data/hg19.chrom.sizes`).

The biological question that will be anwered is "*Is their any potential CNVs can be addressed in the 20 cells by exploring their scRNA-seq data and visualizing in the heatmap?*"

The expected output of this tool will be a heatmap contains the normalized number of reads in each small windows in the whole genome for the 20 cells. There might be some horizontal lighter colored bands which might be potential CNV sites, because these lighter colored horizontal bands means less reads in such small windows.


