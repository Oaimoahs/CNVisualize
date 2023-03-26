from CNVisualize import get_dataframe as m

def test_get_numreads_gen():
    bam = "./../../data/bams/ACGCCAGAGCAGTACG-1_B12_sorted.bam"
    region = "1:500001-1000000"
    expected_numreads = "439"
    numreads = m.get_numreads(bam, region)
    assert numreads == expected_numreads, "The number of reads in %s of this bam is %s not %s" % (region, expected_numreads, numreads)