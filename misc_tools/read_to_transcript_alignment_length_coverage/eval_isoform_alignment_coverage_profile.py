#!/usr/bin/env python3
import argparse
import pysam
import logging
from collections import defaultdict
import time
import math
import statistics


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine BAM read statistics with coverage profiles (from samtools depth)."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument(
        "-d", "--depth", required=True, help="Depth file from 'samtools depth'"
    )
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument(
        "--bins",
        type=int,
        default=100,
        help="Number of bins per transcript for normalized coverage (default: 100)",
    )
    parser.add_argument(
        "--end-fraction",
        type=float,
        default=0.10,
        help="Fraction of bins to use at each end for end-based metrics (default: 0.10)",
    )
    parser.add_argument(
        "--pseudo",
        type=float,
        default=1e-12,
        help="Pseudocount for ratio/contrast stability (default: 1e-12)",
    )
    parser.add_argument(
        "--update-every",
        type=int,
        default=1_000_000,
        help="Log a progress update every N records (default: 1,000,000)",
    )
    parser.add_argument(
        "--log",
        default="INFO",
        help="Logging level: DEBUG, INFO, WARNING, ERROR (default: INFO)",
    )
    return parser.parse_args()


def get_bam_stats_and_lengths(bam_file, update_every):
    """
    Collect read stats per reference and reference lengths from the BAM.
    """
    logging.info(f"Scanning BAM: {bam_file}")
    stats = defaultdict(lambda: {"read_count": 0, "fwd_count": 0, "rc_count": 0})
    nreads = 0
    start = time.time()

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))

        for read in bam.fetch(until_eof=True):
            nreads += 1
            ref = (
                bam.get_reference_name(read.reference_id)
                if read.reference_id >= 0
                else None
            )
            if ref is None:
                continue
            stats[ref]["read_count"] += 1
            qname = read.query_name
            if ":fwd:" in qname:
                stats[ref]["fwd_count"] += 1
            elif ":rc:" in qname:
                stats[ref]["rc_count"] += 1

            if nreads % update_every == 0:
                elapsed = time.time() - start
                logging.info(f"Processed {nreads:,} BAM reads (elapsed {elapsed:.1f}s)")

    # derived metrics
    for ref, d in stats.items():
        rc, fwd, tot = d["rc_count"], d["fwd_count"], d["read_count"]
        d["fwd_rev_ratio"] = (fwd / rc) if rc > 0 else float("inf")
        d["fraction_antisense"] = (rc / tot) if tot > 0 else 0.0

    logging.info(f"Completed BAM scan: {nreads:,} reads, {len(stats)} references")
    return stats, ref_lengths


def compute_end_metrics(norm_bins, end_fraction=0.10, pseudo=1e-12):
    """
    Compute 5' vs 3' end metrics from normalized coverage bins.
    Returns: fiveprime_mean, threeprime_mean, ratio (C5/C3), bias_contrast, log2_ratio
    """
    n = len(norm_bins)
    if n == 0:
        return (math.nan, math.nan, math.nan, math.nan, math.nan)

    # Ensure bins sum to 1 for comparability
    total = sum(norm_bins)
    if total > 0:
        ys = [y / total for y in norm_bins]
    else:
        ys = [0.0] * n

    k = max(1, int(n * end_fraction))
    c5 = sum(ys[:k]) / k
    c3 = sum(ys[-k:]) / k

    ratio = (c5 / c3) if c3 > 0 else (float("inf") if c5 > 0 else 0.0)
    bias_contrast = (c5 - c3) / (c5 + c3 + pseudo)  # bounded [-1,1]
    log2_ratio = math.log2((c5 + pseudo) / (c3 + pseudo))

    return (c5, c3, ratio, bias_contrast, log2_ratio)


def process_depth(
    depth_file, stats, ref_lengths, nbins, end_fraction, pseudo, outfh, update_every
):
    """
    Stream the depth file, bin coverage, and write combined output per transcript.
    """
    logging.info(f"Processing depth file: {depth_file}")
    current_ref = None
    raw_bins = []
    line_counter = 0
    start = time.time()

    def flush(ref, raw_bins):
        if not ref:
            return
        length = ref_lengths.get(ref, 0)
        total = sum(raw_bins)
        norm_bins = [x / total if total > 0 else 0.0 for x in raw_bins]

        # End-based metrics (means, ratio, contrast, log2ratio)
        c5, c3, ratio, contrast, log2ratio = compute_end_metrics(
            norm_bins, end_fraction=end_fraction, pseudo=pseudo
        )

        # Keep older BAM-derived stats
        s = stats.get(
            ref,
            {
                "read_count": 0,
                "fwd_count": 0,
                "rc_count": 0,
                "fwd_rev_ratio": 0,
                "fraction_antisense": 0,
            },
        )

        outfh.write(
            "\t".join(
                [
                    ref,  # transcript_id
                    str(length),  # transcript_length
                    str(s["read_count"]),
                    str(s["fwd_count"]),
                    str(s["rc_count"]),
                    f"{s['fwd_rev_ratio']:.6g}",
                    f"{s['fraction_antisense']:.6g}",
                    f"{ratio:.6g}",  # fiveprime_threeprime_ratio (C5/C3)
                    f"{contrast:.6g}",  # bias_contrast
                    f"{log2ratio:.6g}",  # fiveprime_threeprime_log2ratio
                    f"{c5:.6g}",  # fiveprime_mean
                    f"{c3:.6g}",  # threeprime_mean
                    ",".join(map(str, raw_bins)),
                    ",".join(f"{v:.6g}" for v in norm_bins),
                ]
            )
            + "\n"
        )

    with open(depth_file) as fh:
        for line in fh:
            line_counter += 1
            ref, pos, depth = line.rstrip().split("\t")
            pos, depth = int(pos), int(depth)

            if ref != current_ref:
                flush(current_ref, raw_bins)
                current_ref = ref
                raw_bins = [0] * nbins

            length = ref_lengths.get(ref, 0)
            if length > 0:
                bin_idx = int((pos - 1) / length * nbins)
                bin_idx = min(bin_idx, nbins - 1)
                raw_bins[bin_idx] += depth

            if line_counter % update_every == 0:
                elapsed = time.time() - start
                logging.info(
                    f"Processed {line_counter:,} depth lines (elapsed {elapsed:.1f}s)"
                )

        flush(current_ref, raw_bins)

    logging.info(f"Completed depth processing: {line_counter:,} lines")


def main():
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    stats, ref_lengths = get_bam_stats_and_lengths(args.bam, args.update_every)

    with open(args.output, "w") as out:
        header = [
            "transcript_id",
            "transcript_length",
            "read_count",
            "fwd_count",
            "rc_count",
            "fwd_rev_ratio",
            "fraction_antisense",
            "fiveprime_threeprime_ratio",  # C5/C3
            "bias_contrast",  # (C5 - C3) / (C5 + C3 + pseudo)
            "fiveprime_threeprime_log2ratio",  # log2((C5+p)/(C3+p))
            "fiveprime_mean",
            "threeprime_mean",
            "raw_cov_bins",
            "norm_cov_bins",
        ]
        out.write("\t".join(header) + "\n")

        process_depth(
            args.depth,
            stats,
            ref_lengths,
            args.bins,
            args.end_fraction,
            args.pseudo,
            out,
            args.update_every,
        )

    logging.info(f"Analysis complete. Results written to {args.output}")


if __name__ == "__main__":
    main()
