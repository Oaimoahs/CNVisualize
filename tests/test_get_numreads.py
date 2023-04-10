import CNVisualize.count as m


def test_get_numreads():
    bam = "/Users/shaomiao/Desktop/BIOINF576/CNVisualize/data/sorted_split_bams/ACGCCAGAGCAGTACG-1_B12_sorted.bam"
    region = "1:500001-1000000"
    expected_numreads = "439"
    numreads = m.get_numreads(bam, region)
    assert numreads == expected_numreads, \
        "The number of reads in %s of this bam is %s not %s" % (region, expected_numreads, numreads)
