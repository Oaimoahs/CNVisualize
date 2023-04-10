import CNVisualize.count as m


def test_read_chrom_size():
    path = "/Users/shaomiao/Desktop/BIOINF576/CNVisualize/data/hg19.chrom.sizes"
    expected_chrom_size_dict = {
        "1": 249250621,
        "2": 243199373,
        "3": 198022430,
        "4": 191154276,
        "5": 180915260,
        "6": 171115067,
        "7": 159138663,
        "8": 146364022,
        "9": 141213431,
        "10": 135534747,
        "11": 135006516,
        "12": 133851895,
        "13": 115169878,
        "14": 107349540,
        "15": 102531392,
        "16": 90354753,
        "17": 81195210,
        "18": 78077248,
        "19": 59128983,
        "20": 63025520,
        "21": 48129895,
        "22": 51304566,
        "X": 155270560,
        "Y": 59373566,
    }
    chrom_size_dict = m.read_chrom_size(path)
    assert chrom_size_dict == expected_chrom_size_dict, "The generated dictionary is not the same as the expected one."
