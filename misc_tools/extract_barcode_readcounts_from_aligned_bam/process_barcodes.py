#!/usr/bin/env python3
"""
process_barcodes.py

Process BAM files to extract barcode read counts and match to patient IDs.

Usage:
    python process_barcodes.py <bam_file> <metadata_file> --barcode-tag CB --output-prefix output
"""

import pysam
import pandas as pd
from collections import defaultdict
import re
import sys
import os
import argparse

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def extract_pool_id_from_path(bam_path):
    """Extract pool ID from BAM filename."""
    filename = os.path.basename(bam_path)
    patterns = [
        r'(Neomet_[Pp]ool_?\d+)',
        r'(NeoMet_[Pp]ool_?\d+)',
    ]
    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            pool_id = match.group(1)
            pool_id = re.sub(r'_0+(\d+)', r'_\1', pool_id)
            pool_id = re.sub(r'pool', 'Pool', pool_id, flags=re.IGNORECASE)
            return pool_id
    return "Unknown_Pool"

def normalize_pool_id(pool_id):
    """Normalize pool ID for matching."""
    if not pool_id:
        return None
    pool_id = pool_id.lower()
    pool_id = re.sub(r'_0+(\d+)', r'_\1', pool_id)
    return pool_id

def load_metadata(metadata_file):
    """Load metadata and create barcode to patient mapping."""
    print(f"Loading metadata from {metadata_file}...", file=sys.stderr)
    df = pd.read_csv(metadata_file, sep='\t')
    barcode_map = {}
    pool_patients = defaultdict(set)
    
    for _, row in df.iterrows():
        barcode = row['barcode']
        assignment = row['assignment']
        pool_id = row['pool_id']
        barcode_map[barcode] = (assignment, pool_id)
        pool_patients[pool_id].add(assignment)
    
    print(f"Loaded {len(barcode_map)} barcodes from metadata", file=sys.stderr)
    print(f"Found {len(pool_patients)} pools in metadata", file=sys.stderr)
    
    return barcode_map, pool_patients

def extract_barcodes_from_bam(bam_file, barcode_tag='CB'):
    """Extract barcode read counts from BAM file."""
    print(f"Extracting barcodes from BAM: {bam_file}", file=sys.stderr)
    print(f"Using barcode tag: {barcode_tag}", file=sys.stderr)
    barcode_counts = defaultdict(int)
    
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        total_reads = 0
        reads_with_barcode = 0
        
        # Extract barcodes directly
        for read in bamfile.fetch():
            total_reads += 1
            
            if read.has_tag(barcode_tag):
                barcode = read.get_tag(barcode_tag)
                # Remove suffix like -1
                barcode = barcode.split('-')[0] if isinstance(barcode, str) and '-' in barcode else str(barcode)
                barcode_counts[barcode] += 1
                reads_with_barcode += 1
            
            # Progress indicator for large files
            if total_reads % 1000000 == 0:
                print(f"  Processed {total_reads:,} reads...", file=sys.stderr)
        
        bamfile.close()
        
        print(f"Total reads: {total_reads:,}", file=sys.stderr)
        print(f"Reads with barcode: {reads_with_barcode:,}", file=sys.stderr)
        print(f"Unique barcodes: {len(barcode_counts):,}", file=sys.stderr)
        
    except Exception as e:
        print(f"Error processing BAM: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return {}
    
    return dict(barcode_counts)

def match_barcodes_to_patients(barcode_counts, barcode_map, pool_id, pool_patients):
    """Match barcodes from BAM to patients in metadata."""
    print(f"Matching barcodes to patients for pool: {pool_id}", file=sys.stderr)
    
    # Normalize and find matching pool in metadata
    norm_bam_pool = normalize_pool_id(pool_id)
    matching_pools = []
    
    for meta_pool in pool_patients.keys():
        if normalize_pool_id(meta_pool) == norm_bam_pool:
            matching_pools.append(meta_pool)
    
    if not matching_pools:
        print(f"Warning: Could not match pool '{pool_id}' to metadata", file=sys.stderr)
        print(f"Available pools: {list(pool_patients.keys())}", file=sys.stderr)
        # Try partial matching
        for meta_pool in pool_patients.keys():
            if norm_bam_pool and norm_bam_pool.replace('_', '') in normalize_pool_id(meta_pool).replace('_', ''):
                matching_pools.append(meta_pool)
                print(f"Using partial match: {meta_pool}", file=sys.stderr)
                break
    else:
        print(f"Matched to metadata pool(s): {matching_pools}", file=sys.stderr)
    
    # Process barcodes
    patient_stats = defaultdict(lambda: {'barcode_count': 0, 'total_reads': 0})
    barcode_details = []
    unmatched =
