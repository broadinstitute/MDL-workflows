#author: cgeorges@broadinstitute.org

#!/usr/bin/env python3

import argparse
import pickle
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple
import pysam


def list_contigs(bam_path: str) -> List[Tuple[str, int]]:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return list(zip(bam.header.references, bam.header.lengths))


def count_lengths_contig(bam_path: str, contig: str) -> Counter:
    c = Counter()
    # threads=1 to avoid oversubscription across processes
    with pysam.AlignmentFile(bam_path, "rb", threads=1) as bam:
        for r in bam.fetch(contig):
            ql = r.query_length
            if ql:
                c[ql] += 1
    return c


def count_lengths_single_pass(bam_path: str, bgzf_threads: int) -> Counter:
    c = Counter()
    with pysam.AlignmentFile(bam_path, "rb", threads=max(1, bgzf_threads)) as bam:
        for r in bam.fetch(until_eof=True):
            ql = r.query_length
            if ql:
                c[ql] += 1
    return c


def run_parallel_mapped(bam_path: str, workers: int) -> Counter:
    contigs = list_contigs(bam_path)
    if not contigs:
        raise RuntimeError("No contigs in header.")
    # largest-first; executor handles dynamic scheduling
    ctg_names = [n for n, L in sorted(contigs, key=lambda x: x[1], reverse=True)]
    total = Counter()
    with ProcessPoolExecutor(max_workers=max(1, workers)) as ex:
        futs = [ex.submit(count_lengths_contig, bam_path, ctg) for ctg in ctg_names]
        for f in as_completed(futs):
            total.update(f.result())
    return total


def main():
    ap = argparse.ArgumentParser(
        description="Read-length histogram from BAM/CRAM. "
                    "Parallel per-contig for mapped-only. Single-pass when including unmapped."
    )
    ap.add_argument("-i", "--bam", required=True, help="Input BAM/CRAM")
    ap.add_argument("-o", "--out", required=True, help="Output pickle path")
    ap.add_argument("--include-unmapped", action="store_true",
                    help="Single-pass over entire file (counts mapped and unmapped).")
    ap.add_argument("-w", "--workers", type=int, default=4,
                    help="Processes for per-contig mode (ignored if --include-unmapped).")
    ap.add_argument("-t", "--bgzf-threads", type=int, default=4,
                    help="BGZF decompression threads for single-pass mode (only used with --include-unmapped).")
    args = ap.parse_args()

    if args.include_unmapped:
        counts = count_lengths_single_pass(args.bam, args.bgzf_threads)
    else:
        counts = run_parallel_mapped(args.bam, args.workers)

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with outp.open("wb") as fh:
        pickle.dump(counts, fh, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()


