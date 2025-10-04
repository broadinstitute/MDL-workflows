#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict

def count_per_cell(bam_file, cb_tag="CB", umi_tag="XM", threads=1):
    bam = pysam.AlignmentFile(bam_file, "rb")
    read_counts = defaultdict(int)
    umi_counts = defaultdict(set)

    for read in bam.fetch(until_eof=True):
        if read.has_tag(cb_tag):
            cb = read.get_tag(cb_tag)
            read_counts[cb] += 1
            if read.has_tag(umi_tag):
                umi = read.get_tag(umi_tag)
                umi_counts[cb].add(umi)

    bam.close()
    umi_final = {cb: len(umis) for cb, umis in umi_counts.items()}
    return read_counts, umi_final


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads and UMIs per cell barcode from BAM")
    parser.add_argument("bam", help="Input BAM file (must be indexed)")
    parser.add_argument("--cb_tag", default="CB", help="Cell barcode tag (default: CB)")
    parser.add_argument("--umi_tag", default="XM", help="UMI tag (default: XM)")
    parser.add_argument("--sample_id", required=True, help="Sample identifier to label results")
    parser.add_argument("-o", "--output", default="counts.tsv", help="Output TSV file")

    args = parser.parse_args()

    read_counts, umi_counts = count_per_cell(args.bam, args.cb_tag, args.umi_tag)

    with open(args.output, "w") as out:
        out.write("sample_id\tcell_barcode\treads\tumis\n")
        for cb in sorted(read_counts.keys()):
            reads = read_counts[cb]
            umis = umi_counts.get(cb, 0)
            out.write(f"{args.sample_id}\t{cb}\t{reads}\t{umis}\n")
