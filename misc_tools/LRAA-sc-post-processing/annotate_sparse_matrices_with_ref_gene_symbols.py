#!/usr/bin/env python3

'''
Usage:
python3 annotate_sparse_matrices_with_ref_gene_symbols.py \
    --reference_gtf reference.gtf \
    --gene_sparse_dir gene_matrix_dir \
    --isoform_sparse_dir isoform_matrix_dir \
    --output_gene_dir gene_matrix_annotated \
    --output_isoform_dir isoform_matrix_annotated \
    --gene_mappings_output gene_symbol_mappings.tsv
'''

import sys
import os
import re
import gzip
import argparse
import pandas as pd
import shutil

def open_file(filename):
    """Open file handling both regular and gzipped files"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def extract_gene_symbols_from_gtf(gtf_file):
    """Extract gene symbols from reference GTF file"""
    gene_id_to_symbol = {}
    
    print(f"Processing reference GTF for gene symbols: {gtf_file}")
    
    try:
        with open_file(gtf_file) as fh:
            for line_num, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                    
                feature_type = fields[2]
                if feature_type not in ['gene', 'transcript']:
                    continue
                    
                attributes = fields[8]
                
                gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
                gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)
                transcript_id_match = re.search(r'transcript_id "([^"]+)"', attributes)
                
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                    gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                    gene_id_to_symbol[gene_id] = gene_name
                    
                    # Also map transcript ID to gene symbol for transcript-level matrices
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)
                        gene_id_to_symbol[transcript_id] = gene_name
                
                # Progress indicator for large files
                if line_num % 100000 == 0:
                    print(f"Processed {line_num} lines, found {len(gene_id_to_symbol)} mappings")
                    
    except Exception as e:
        print(f"Error processing GTF file: {e}")
        return {}

    print(f"Extracted {len(gene_id_to_symbol)} gene/transcript symbol mappings")
    return gene_id_to_symbol

def annotate_sparse_matrix(matrix_dir, gene_mapping, output_dir):
    """Annotate sparse matrix features with gene symbols"""
    
    print(f"Processing matrix directory: {matrix_dir}")
    
    if not os.path.exists(matrix_dir):
        print(f"Error: Matrix directory {matrix_dir} does not exist")
        return False
    
    # Find features file
    features_file = None
    for filename in ['features.tsv', 'genes.tsv', 'features.tsv.gz', 'genes.tsv.gz']:
        potential_file = os.path.join(matrix_dir, filename)
        if os.path.exists(potential_file):
            features_file = potential_file
            break
    
    if features_file is None:
        print(f"Error: No features file found in {matrix_dir}")
        print(f"Available files: {os.listdir(matrix_dir) if os.path.exists(matrix_dir) else 'Directory does not exist'}")
        return False
    
    print(f"Reading features from: {features_file}")
    
    # Read features
    try:
        if features_file.endswith('.gz'):
            features = pd.read_csv(features_file, sep='\t', header=None, compression='gzip')
        else:
            features = pd.read_csv(features_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading features file: {e}")
        return False
    
    print(f"Original features shape: {features.shape}")
    print(f"First few features:\n{features.head()}")
    
    # Ensure we have at least 2 columns
    if len(features.columns) == 1:
        features[1] = features[0]  # Duplicate first column as feature name
    elif len(features.columns) == 0:
        print("Error: Features file is empty")
        return False
    
    # Set column names
    features.columns = ['feature_id', 'feature_name'] + [f'col_{i}' for i in range(2, len(features.columns))]
    
    # Annotate with gene symbols
    annotated_features = features.copy()
    
    def get_gene_symbol(feature_id):
        """Get gene symbol for a feature ID"""
        if pd.isna(feature_id) or feature_id == '':
            return feature_id
            
        feature_id_str = str(feature_id)
        
        if feature_id_str in gene_mapping:
            gene_symbol = gene_mapping[feature_id_str]
            return f"{gene_symbol}^{feature_id_str}"
        else:
            # For transcript IDs, try to extract gene ID
            if feature_id_str.startswith('ENST'):
                # This is a transcript ID, we might not have direct mapping
                # Keep original for now
                return feature_id_str
            return feature_id_str
    
    annotated_features['annotated_name'] = annotated_features['feature_id'].apply(get_gene_symbol)
    
    # Create output directory
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Error creating output directory {output_dir}: {e}")
        return False
    
    # Copy all files from original directory except features files
    try:
        for file in os.listdir(matrix_dir):
            if file not in ['features.tsv', 'genes.tsv', 'features.tsv.gz', 'genes.tsv.gz']:
                src_file = os.path.join(matrix_dir, file)
                dst_file = os.path.join(output_dir, file)
                if os.path.isfile(src_file):
                    shutil.copy2(src_file, dst_file)
                elif os.path.isdir(src_file):
                    shutil.copytree(src_file, dst_file, dirs_exist_ok=True)
    except Exception as e:
        print(f"Warning: Error copying files: {e}")
    
    # Write annotated features
    output_features_file = os.path.join(output_dir, 'features.tsv')
    output_df = pd.DataFrame({
        'feature_id': annotated_features['feature_id'],
        'feature_name': annotated_features['annotated_name']
    })
    
    try:
        output_df.to_csv(output_features_file, sep='\t', header=False, index=False)
        print(f"Annotated features saved to: {output_features_file}")
        print(f"Annotated features shape: {output_df.shape}")
        print(f"Sample annotated features:\n{output_df.head()}")
        return True
    except Exception as e:
        print(f"Error writing annotated features: {e}")
        return False

def save_gene_symbol_mappings(gene_mapping, output_file):
    """Save gene symbol mappings to file"""
    try:
        with open(output_file, "w") as out_fh:
            out_fh.write("feature_id\tgene_symbol\n")
            for feature_id, symbol in gene_mapping.items():
                out_fh.write(f"{feature_id}\t{symbol}\n")
        print(f"Gene symbol mapping saved to {output_file}")
        return True
    except Exception as e:
        print(f"Error saving gene mappings: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Annotate sparse matrices with gene symbols from reference GTF",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "--reference_gtf", 
        required=True,
        help="Reference GTF file (e.g., GENCODE)"
    )
    
    parser.add_argument(
        "--gene_sparse_dir",
        required=True, 
        help="Gene sparse matrix directory"
    )
    
    parser.add_argument(
        "--isoform_sparse_dir",
        required=True,
        help="Isoform sparse matrix directory" 
    )
    
    parser.add_argument(
        "--output_gene_dir",
        required=True,
        help="Output directory for annotated gene sparse matrix"
    )
    
    parser.add_argument(
        "--output_isoform_dir", 
        required=True,
        help="Output directory for annotated isoform sparse matrix"
    )
    
    parser.add_argument(
        "--gene_mappings_output",
        default="gene_symbol_mappings.tsv",
        help="Output file for gene symbol mappings (default: gene_symbol_mappings.tsv)"
    )
    
    args = parser.parse_args()
    
    # Validate input files and directories
    if not os.path.exists(args.reference_gtf):
        print(f"Error: Reference GTF file does not exist: {args.reference_gtf}")
        sys.exit(1)
    
    if not os.path.exists(args.gene_sparse_dir):
        print(f"Error: Gene sparse directory does not exist: {args.gene_sparse_dir}")
        sys.exit(1)
        
    if not os.path.exists(args.isoform_sparse_dir):
        print(f"Error: Isoform sparse directory does not exist: {args.isoform_sparse_dir}")
        sys.exit(1)
    
    # Extract gene symbols from reference GTF
    print("=" * 60)
    print("EXTRACTING GENE SYMBOLS FROM REFERENCE GTF")
    print("=" * 60)
    
    gene_id_to_symbol = extract_gene_symbols_from_gtf(args.reference_gtf)
    
    if not gene_id_to_symbol:
        print("Error: No gene symbols extracted from GTF")
        sys.exit(1)
    
    # Save gene mappings
    if not save_gene_symbol_mappings(gene_id_to_symbol, args.gene_mappings_output):
        print("Error: Failed to save gene mappings")
        sys.exit(1)
    
    # Process gene sparse matrix
    print("\n" + "=" * 60)
    print("ANNOTATING GENE SPARSE MATRIX")
    print("=" * 60)
    
    success_gene = annotate_sparse_matrix(
        args.gene_sparse_dir, 
        gene_id_to_symbol, 
        args.output_gene_dir
    )
    
    # Process isoform sparse matrix  
    print("\n" + "=" * 60)
    print("ANNOTATING ISOFORM SPARSE MATRIX")
    print("=" * 60)
    
    success_isoform = annotate_sparse_matrix(
        args.isoform_sparse_dir,
        gene_id_to_symbol,
        args.output_isoform_dir
    )
    
    # Summary
    print("\n" + "=" * 60)
    print("ANNOTATION SUMMARY")
    print("=" * 60)
    print(f"Gene symbols extracted: {len(gene_id_to_symbol)}")
    print(f"Gene matrix annotation: {'SUCCESS' if success_gene else 'FAILED'}")
    print(f"Isoform matrix annotation: {'SUCCESS' if success_isoform else 'FAILED'}")
    
    if success_gene and success_isoform:
        print("\nAll annotations completed successfully!")
        sys.exit(0)
    else:
        print("\nSome annotations failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()
