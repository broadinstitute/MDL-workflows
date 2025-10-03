#!/usr/bin/env python3
"""
Merge barcode counting results from multiple pool BAM files.
"""

import pandas as pd
import sys
import argparse
from collections import defaultdict
import os

def merge_patient_stats(patient_stats_files):
    """Merge patient statistics from multiple pools."""
    print("Merging patient stats...", file=sys.stderr)
    all_patient_stats = defaultdict(lambda: {'barcode_count': 0, 'total_reads': 0, 'pools': []})
    
    for file in patient_stats_files:
        if not file or not os.path.exists(file):
            continue
            
        df = pd.read_csv(file, sep='\t')
        pool_name = os.path.basename(file).replace('_patient_stats.tsv', '')
        
        for _, row in df.iterrows():
            patient_id = row['patient_id']
            all_patient_stats[patient_id]['barcode_count'] += int(row['barcode_count'])
            all_patient_stats[patient_id]['total_reads'] += int(row['total_reads'])
            all_patient_stats[patient_id]['pools'].append(pool_name)
    
    # Create merged patient summary
    summary_data = []
    for patient_id, stats in sorted(all_patient_stats.items()):
        summary_data.append({
            'patient_id': patient_id,
            'barcode_count': stats['barcode_count'],
            'total_reads': stats['total_reads'],
            'avg_reads_per_barcode': round(stats['total_reads'] / stats['barcode_count'], 2) if stats['barcode_count'] > 0 else 0,
            'num_pools': len(set(stats['pools'])),
            'pools': ','.join(sorted(set(stats['pools'])))
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('total_reads', ascending=False)
    
    print(f"Total patients: {len(summary_df)}", file=sys.stderr)
    print(f"Total barcodes: {summary_df['barcode_count'].sum():,}", file=sys.stderr)
    print(f"Total reads: {summary_df['total_reads'].sum():,}", file=sys.stderr)
    
    return summary_df

def merge_barcode_details(barcode_files):
    """Combine all barcode details from multiple pools."""
    print("\nCombining barcode details...", file=sys.stderr)
    barcode_dfs = []
    
    for file in barcode_files:
        if not file or not os.path.exists(file):
            continue
            
        df = pd.read_csv(file, sep='\t')
        if not df.empty:
            barcode_dfs.append(df)
    
    if barcode_dfs:
        all_barcodes_df = pd.concat(barcode_dfs, ignore_index=True)
        all_barcodes_df = all_barcodes_df.sort_values(['patient_id', 'read_count'], ascending=[True, False])
        print(f"Total barcode entries: {len(all_barcodes_df):,}", file=sys.stderr)
        return all_barcodes_df
    else:
        print("No barcode details found", file=sys.stderr)
        return pd.DataFrame(columns=['barcode', 'patient_id', 'read_count', 'match_type', 'pool_id'])

def merge_unmatched_barcodes(unmatched_files):
    """Combine all unmatched barcodes from multiple pools."""
    print("\nCombining unmatched barcodes...", file=sys.stderr)
    unmatched_dfs = []
    
    for file in unmatched_files:
        if not file or not os.path.exists(file):
            continue
            
        df = pd.read_csv(file, sep='\t')
        if not df.empty:
            pool_name = os.path.basename(file).replace('_unmatched.tsv', '')
            df['pool'] = pool_name
            unmatched_dfs.append(df)
    
    if unmatched_dfs:
        all_unmatched_df = pd.concat(unmatched_dfs, ignore_index=True)
        all_unmatched_df = all_unmatched_df.sort_values('read_count', ascending=False)
        print(f"Total unmatched barcodes: {len(all_unmatched_df):,}", file=sys.stderr)
        print(f"Total unmatched reads: {all_unmatched_df['read_count'].sum():,}", file=sys.stderr)
        return all_unmatched_df
    else:
        print("No unmatched barcodes found", file=sys.stderr)
        return pd.DataFrame(columns=['barcode', 'read_count', 'pool'])

def main():
    parser = argparse.ArgumentParser(description='Merge barcode counting results from multiple pools')
    parser.add_argument('--patient-stats', nargs='+', required=True, 
                       help='Patient stats TSV files from each pool')
    parser.add_argument('--barcode-details', nargs='+', required=True,
                       help='Barcode details TSV files from each pool')
    parser.add_argument('--unmatched', nargs='+', required=True,
                       help='Unmatched barcodes TSV files from each pool')
    parser.add_argument('--output-prefix', default='merged',
                       help='Output file prefix (default: merged)')
    
    args = parser.parse_args()
    
    print(f"{'='*60}", file=sys.stderr)
    print(f"MERGING RESULTS", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    print(f"Patient stats files: {len(args.patient_stats)}", file=sys.stderr)
    print(f"Barcode details files: {len(args.barcode_details)}", file=sys.stderr)
    print(f"Unmatched files: {len(args.unmatched)}", file=sys.stderr)
    
    # Merge patient stats
    patient_summary_df = merge_patient_stats(args.patient_stats)
    patient_output = f"{args.output_prefix}_patient_summary.tsv"
    patient_summary_df.to_csv(patient_output, sep='\t', index=False)
    print(f"\nSaved: {patient_output}", file=sys.stderr)
    
    # Merge barcode details
    all_barcodes_df = merge_barcode_details(args.barcode_details)
    barcodes_output = f"{args.output_prefix}_all_barcodes.tsv"
    all_barcodes_df.to_csv(barcodes_output, sep='\t', index=False)
    print(f"Saved: {barcodes_output}", file=sys.stderr)
    
    # Merge unmatched barcodes
    unmatched_df = merge_unmatched_barcodes(args.unmatched)
    unmatched_output = f"{args.output_prefix}_unmatched.tsv"
    unmatched_df.to_csv(unmatched_output, sep='\t', index=False)
    print(f"Saved: {unmatched_output}", file=sys.stderr)
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"MERGE COMPLETE", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    # Print top patients
    if not patient_summary_df.empty:
        print(f"\nTop 10 patients by read count:", file=sys.stderr)
        top10 = patient_summary_df.head(10)[['patient_id', 'barcode_count', 'total_reads', 'num_pools']]
        print(top10.to_string(index=False), file=sys.stderr)

if __name__ == "__main__":
    main()
