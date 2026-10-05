#!/usr/bin/env python3

import argparse
import multiprocessing as mp
from array import array
from collections import defaultdict

import numpy as np
import pysam

# UMIs are packed as base-5 integers (N=4) next to a barcode id, 8 bytes per read instead of
# a Python str in a set (~100 bytes). UMIs longer than MAX_PACKED_UMI go to a str-set fallback.
UMI_TABLE = str.maketrans("ACGTN", "01234")
MAX_PACKED_UMI = 13  # 5**13 < 2**32


def count_contigs(args):
    """Count reads/UMIs for the records on `contigs` (None = whole file, via until_eof)."""
    bam_file, contigs, cb_tag, umi_tag, threads, count_umis, max_reads = args
    bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)
    read_counts = defaultdict(int)
    cb_ids = {}
    packed = array("Q")
    fallback = defaultdict(set)

    if contigs is None:
        streams = [bam.fetch(until_eof=True)]
    elif contigs == ["*"]:
        # records with no coordinate (unplaced unmapped reads), which have no contig
        streams = [bam.fetch("*")]
    else:
        streams = [bam.fetch(contig=c) for c in contigs]

    n = 0
    for stream in streams:
        for read in stream:
            n += 1
            if max_reads and n > max_reads:
                break
            try:
                cb = read.get_tag(cb_tag)
            except KeyError:
                continue
            read_counts[cb] += 1
            if count_umis:
                try:
                    umi = read.get_tag(umi_tag)
                except KeyError:
                    continue
                cid = cb_ids.get(cb)
                if cid is None:
                    cid = cb_ids[cb] = len(cb_ids)
                if len(umi) <= MAX_PACKED_UMI:
                    try:
                        packed.append((cid << 32) | int(umi.translate(UMI_TABLE), 5))
                        continue
                    except ValueError:
                        pass
                fallback[cb].add(umi)
    bam.close()
    unique = np.unique(np.frombuffer(packed, dtype=np.uint64)) if len(packed) else np.empty(0, np.uint64)
    return dict(read_counts), list(cb_ids), unique, dict(fallback)



def plan_shards(bam_file, n_shards):
    """Greedy-balance contigs by record count (from the .bai) into n_shards lists."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    stats = {s.contig: s.mapped + s.unmapped for s in bam.get_index_statistics()}
    nocoord = bam.nocoordinate
    bam.close()
    bins = [[0, []] for _ in range(n_shards)]
    for contig, size in sorted(stats.items(), key=lambda kv: -kv[1]):
        if size == 0:
            continue
        smallest = min(bins, key=lambda b: b[0])
        smallest[0] += size
        smallest[1].append(contig)
    shards = [b[1] for b in bins if b[1]]
    if nocoord:
        shards.append(["*"])
    return shards


def count_per_cell(bam_file, cb_tag="CB", umi_tag="XM", threads=1, processes=1,
                   count_umis=True, max_reads=0):
    if processes <= 1 or max_reads:
        jobs = [(bam_file, None, cb_tag, umi_tag, threads, count_umis, max_reads)]
    else:
        jobs = [(bam_file, s, cb_tag, umi_tag, threads, count_umis, 0)
                for s in plan_shards(bam_file, processes)]

    if len(jobs) == 1:
        results = [count_contigs(jobs[0])]
    else:
        with mp.Pool(min(processes, len(jobs))) as pool:
            results = pool.map(count_contigs, jobs, chunksize=1)

    read_counts = defaultdict(int)
    global_ids = {}
    chunks = []
    fallback = defaultdict(set)
    for rc, cb_list, unique, fb in results:
        for cb, n in rc.items():
            read_counts[cb] += n
        # remap worker-local barcode ids to global ids so identical (barcode, UMI) pairs collapse
        remap = np.array([global_ids.setdefault(cb, len(global_ids)) for cb in cb_list], dtype=np.uint64)
        if len(unique):
            chunks.append((remap[unique >> np.uint64(32)] << np.uint64(32)) | (unique & np.uint64(0xFFFFFFFF)))
        for cb, umis in fb.items():
            fallback[cb] |= umis
    umi_final = defaultdict(int)
    if chunks:
        merged = np.unique(np.concatenate(chunks))
        per_cb = np.bincount((merged >> np.uint64(32)).astype(np.int64), minlength=len(global_ids))
        for cb, cid in global_ids.items():
            umi_final[cb] = int(per_cb[cid])
    for cb, umis in fallback.items():
        umi_final[cb] += len(umis)
    return read_counts, dict(umi_final)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count reads and UMIs per cell barcode from BAM")
    parser.add_argument("bam", help="Input BAM file (indexed if --processes > 1)")
    parser.add_argument("--sample_id", required=True, help="Sample identifier to label results")
    parser.add_argument("--cb_tag", default="CB", help="Cell barcode tag (default: CB)")
    parser.add_argument("--umi_tag", default="XM", help="UMI tag (default: XM)")
    parser.add_argument("--processes", type=int, default=1,
                        help="Parallel workers, each handling a share of the contigs (needs .bai)")
    parser.add_argument("--threads", type=int, default=1, help="BGZF decompression threads per worker")
    parser.add_argument("--skip_umi", action="store_true",
                        help="Do not count UMIs (much lower memory); umis column is 0")
    parser.add_argument("--max_reads", type=int, default=0, help="Stop after N records (testing only)")
    parser.add_argument("-o", "--output", default="counts.tsv", help="Output TSV file")

    args = parser.parse_args()

    read_counts, umi_counts = count_per_cell(
        args.bam, args.cb_tag, args.umi_tag, threads=args.threads, processes=args.processes,
        count_umis=not args.skip_umi, max_reads=args.max_reads,
    )

    with open(args.output, "w") as out:
        out.write("sample_id\tcell_barcode\treads\tumis\n")
        for cb in sorted(read_counts.keys()):
            out.write(f"{args.sample_id}\t{cb}\t{read_counts[cb]}\t{umi_counts.get(cb, 0)}\n")
