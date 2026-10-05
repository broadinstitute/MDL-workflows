#!/usr/bin/env python3
"""Sum per-barcode read counts across several process_barcodes.py outputs.

Writes a barcode/post_count TSV, the default input format of assign_barcode_groups.py.
"""

import argparse
import csv
from collections import defaultdict

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("counts", nargs="+", help="TSVs from process_barcodes.py")
parser.add_argument("-o", "--output", default="merged_counts.tsv")
args = parser.parse_args()

totals = defaultdict(int)
for path in args.counts:
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            totals[row["cell_barcode"]] += int(row["reads"])

with open(args.output, "w") as out:
    out.write("barcode\tpost_count\n")
    for barcode, count in sorted(totals.items(), key=lambda t: -t[1]):
        out.write(f"{barcode}\t{count}\n")
