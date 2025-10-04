#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def process_chunk(bam_file, start, end, cb_tag="CB", umi_tag="XM"):
    bam = pysam.AlignmentFile(bam_file, "rb")
    read_counts = defaultdict(int)
    umi_counts = defaultdict(set)

    for read in bam.fetch(until_eof=True, start=start, end=end):
        if read.has_tag(cb_tag):
            cb = read.get_tag(cb_tag)
            read_counts[cb] += 1
            if read.has_tag(umi_tag):
                umi = read.get_tag(umi_tag)
                umi_counts[cb].add(umi)

    bam.close()
    return read_counts, umi_counts

def merge_counts(results):
    final_reads = defaultdict(int)
    final_umis = defaultdict(set)
    for r_counts, u_counts in results:
        for cb, cnt in r_counts.items():
            final_reads[cb] += cnt
        for cb, umis in u_counts.items():
            final_umis[cb].update(umis)
    return final_reads, {cb: len(umis) for cb, umis in final_umis.items()}

def count_per_cell(bam_file, cb_tag="CB", umi_tag="XM", threads=4):
    bam = pysam.AlignmentFile(bam_file, "rb")
    references = bam.references
    lengths = bam.lengths
    bam.close()

    tasks = [(bam_file, 0, length, cb_tag, umi_tag) for length in lengths]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = list(executor.map(lambda args: process_chunk(*args), tasks))

    return merge_counts(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads and UMIs per cell barcode from BAM")
    parser.add_argument("bam", help="Input BAM file (must be indexed)")
    parser.add_argument("--cb_tag", default="CB", help="Cell barcode tag (default: CB)")
    parser.add_argument("--umi_tag", default="XM", help="UMI tag (default: XM)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("-o", "--output", default="cell_counts.tsv", help="Output TSV file")

    args = parser.parse_args()
    read_counts, umi_counts = count_per_cell(args.bam, args.cb_tag, args.umi_tag, args.threads)

    with open(args.output, "w") as out:
        out.write("cell_barcode\treads\tumis\n")
        for cb in sorted(read_counts.keys()):
            reads = read_counts[cb]
            umis = umi_counts.get(cb, 0)
            out.write(f"{cb}\t{reads}\t{umis}\n")
