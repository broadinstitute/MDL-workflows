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
    unmatched = []
    
    for barcode, read_count in barcode_counts.items():
        matched = False
        patient_id = None
        match_type = "direct"
        
        # Try direct match
        if barcode in barcode_map:
            patient_id, meta_pool = barcode_map[barcode]
            if meta_pool in matching_pools or not matching_pools:
                matched = True
        
        # Try reverse complement
        if not matched:
            rev_comp = reverse_complement(barcode)
            if rev_comp in barcode_map:
                patient_id, meta_pool = barcode_map[rev_comp]
                if meta_pool in matching_pools or not matching_pools:
                    matched = True
                    match_type = "reverse_complement"
        
        if matched and patient_id:
            patient_stats[patient_id]['barcode_count'] += 1
            patient_stats[patient_id]['total_reads'] += read_count
            barcode_details.append({
                'barcode': barcode,
                'patient_id': patient_id,
                'read_count': read_count,
                'match_type': match_type,
                'pool_id': pool_id
            })
        else:
            unmatched.append({
                'barcode': barcode,
                'read_count': read_count
            })
    
    return dict(patient_stats), barcode_details, unmatched

def process_bam(bam_file, metadata_file, barcode_tag, output_prefix):
    """Main processing function."""
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"Processing BAM: {os.path.basename(bam_file)}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    # Extract pool ID
    pool_id = extract_pool_id_from_path(bam_file)
    print(f"Detected pool ID: {pool_id}", file=sys.stderr)
    
    # Load metadata
    barcode_map, pool_patients = load_metadata(metadata_file)
    
    # Extract barcodes from BAM
    barcode_counts = extract_barcodes_from_bam(bam_file, barcode_tag)
    
    if not barcode_counts:
        print("No barcodes found! Creating empty outputs.", file=sys.stderr)
        # Create empty outputs
        pd.DataFrame(columns=['patient_id', 'barcode_count', 'total_reads', 'avg_reads_per_barcode']).to_csv(
            f"{output_prefix}_patient_stats.tsv", sep='\t', index=False)
        pd.DataFrame(columns=['barcode', 'patient_id', 'read_count', 'match_type', 'pool_id']).to_csv(
            f"{output_prefix}_barcode_details.tsv", sep='\t', index=False)
        pd.DataFrame(columns=['barcode', 'read_count']).to_csv(
            f"{output_prefix}_unmatched.tsv", sep='\t', index=False)
        return
    
    # Match barcodes to patients
    patient_stats, barcode_details, unmatched = match_barcodes_to_patients(
        barcode_counts, barcode_map, pool_id, pool_patients
    )
    
    # Save patient stats (aggregated per patient)
    patient_df = pd.DataFrame([
        {
            'patient_id': pid,
            'total_umis': stats['barcode_count'],  # Each unique barcode = 1 UMI
            'total_reads': stats['total_reads'],
            'avg_reads_per_umi': round(stats['total_reads'] / stats['barcode_count'], 2) if stats['barcode_count'] > 0 else 0
        }
        for pid, stats in sorted(patient_stats.items())
    ])
    patient_df.to_csv(f"{output_prefix}_patient_stats.tsv", sep='\t', index=False)
    print(f"Saved patient stats: {output_prefix}_patient_stats.tsv", file=sys.stderr)
    
    # Save barcode details (one row per barcode with patient assignment)
    barcode_df = pd.DataFrame(barcode_details)
    if not barcode_df.empty:
        barcode_df = barcode_df.sort_values(['patient_id', 'read_count'], ascending=[True, False])
    barcode_df.to_csv(f"{output_prefix}_barcode_details.tsv", sep='\t', index=False)
    print(f"Saved barcode details: {output_prefix}_barcode_details.tsv", file=sys.stderr)
    
    # Save unmatched barcodes
    unmatched_df = pd.DataFrame(unmatched)
    if not unmatched_df.empty:
        unmatched_df = unmatched_df.sort_values('read_count', ascending=False)
    unmatched_df.to_csv(f"{output_prefix}_unmatched.tsv", sep='\t', index=False)
    print(f"Saved unmatched barcodes: {output_prefix}_unmatched.tsv", file=sys.stderr)
    
    # Print summary
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"SUMMARY", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Pool: {pool_id}", file=sys.stderr)
    print(f"Patients found: {len(patient_stats)}", file=sys.stderr)
    print(f"Barcodes matched: {len(barcode_details):,}", file=sys.stderr)
    print(f"Total reads matched: {sum(b['read_count'] for b in barcode_details):,}", file=sys.stderr)
    print(f"Barcodes unmatched: {len(unmatched):,}", file=sys.stderr)
    print(f"Total reads unmatched: {sum(u['read_count'] for u in unmatched):,}", file=sys.stderr)
    
    if patient_stats:
        print(f"\nPer-patient breakdown:", file=sys.stderr)
        for pid, stats in sorted(patient_stats.items(), key=lambda x: x[1]['total_reads'], reverse=True):
            print(f"  {pid}: {stats['barcode_count']} barcodes, {stats['total_reads']:,} reads", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Process BAM file to extract barcode counts per patient')
    parser.add_argument('bam_file', help='Input BAM file')
    parser.add_argument('metadata_file', help='Metadata TSV file')
    parser.add_argument('--barcode-tag', default='CB', help='BAM tag containing barcodes (default: CB)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    process_bam(args.bam_file, args.metadata_file, args.barcode_tag, args.output_prefix)
    print("\nDone!", file=sys.stderr)

if __name__ == "__main__":
    main()
