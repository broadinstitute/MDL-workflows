import pysam
import pandas as pd
from pathlib import Path
import sys
import gc
import os

def get_required_read_ids(fusion_tsv_path):
    """
    Extract all unique read IDs from the fusion TSV file.
    This way we only need to store tags for reads we actually need.
    """
    print(f"Reading fusion TSV to get required read IDs: {fusion_tsv_path}")
    df = pd.read_csv(fusion_tsv_path, sep="\t")
    print(f"Fusion TSV shape: {df.shape}")
    
    # Get the evidence column (assuming it's the 10th column, index 9)
    evidence_col = df.columns[9]
    print(f"Evidence column: {evidence_col}")
    
    # Extract all read IDs from all rows
    all_read_ids = set()
    for idx, row in df.iterrows():
        if pd.notna(row[evidence_col]):
            reads = str(row[evidence_col]).split(",")
            for read_id in reads:
                read_id = read_id.strip()
                if read_id:
                    all_read_ids.add(read_id)
    
    print(f"Found {len(all_read_ids)} unique read IDs needed from BAM file")
    return all_read_ids

def extract_tags_from_bam_selective(bam_path, required_read_ids):
    """
    Only extract tags for reads that we actually need.
    This dramatically reduces memory usage.
    """
    print(f"Processing BAM file: {bam_path}")
    print(f"Looking for {len(required_read_ids)} specific reads")
    
    tag_dict = {}
    reads_processed = 0
    reads_found = 0
    
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                reads_processed += 1
                
                # Only process reads we actually need
                if read.query_name in required_read_ids:
                    xm = read.get_tag("XM") if read.has_tag("XM") else None
                    cb = read.get_tag("CB") if read.has_tag("CB") else None
                    tag_dict[read.query_name] = (xm, cb)
                    reads_found += 1
                
                # Progress reporting
                if reads_processed % 500000 == 0:
                    print(f"Processed {reads_processed:,} reads, found {reads_found:,} target reads...")
                    # Force garbage collection
                    gc.collect()
                
                # Early exit if we found all required reads
                if reads_found == len(required_read_ids):
                    print(f"Found all required reads after processing {reads_processed:,} total reads")
                    break
                    
    except Exception as e:
        print(f"Error processing BAM file: {e}")
        raise
    
    print(f"Total reads processed: {reads_processed:,}")
    print(f"Target reads found: {reads_found:,} out of {len(required_read_ids):,} required")
    
    # Report missing reads
    missing_reads = required_read_ids - set(tag_dict.keys())
    if missing_reads:
        print(f"Warning: {len(missing_reads)} reads not found in BAM file")
        if len(missing_reads) <= 10:
            print(f"Missing reads: {list(missing_reads)}")
    
    return tag_dict

def update_fusion_tsv(tsv_path, out_path, tag_dict):
    print(f"Reading fusion TSV: {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"Fusion TSV shape: {df.shape}")

    # assuming evidence reads are in the 10th column (index 9)
    evidence_col = df.columns[9]
    print(f"Evidence column: {evidence_col}")

    def get_cb_list(reads_str):
        if pd.isna(reads_str):
            return ""
        reads = str(reads_str).split(",")
        cb_list = []
        for r in reads:
            r = r.strip()
            if r in tag_dict:
                cb = tag_dict[r][1]
                cb_list.append(str(cb) if cb is not None else "None")
            else:
                cb_list.append("None")
        return ",".join(cb_list)

    def get_xm_list(reads_str):
        if pd.isna(reads_str):
            return ""
        reads = str(reads_str).split(",")
        xm_list = []
        for r in reads:
            r = r.strip()
            if r in tag_dict:
                xm = tag_dict[r][0]
                xm_list.append(str(xm) if xm is not None else "None")
            else:
                xm_list.append("None")
        return ",".join(xm_list)

    print("Adding CB_list column...")
    df["CB_list"] = df[evidence_col].apply(get_cb_list)
    
    print("Adding XM_list column...")
    df["XM_list"] = df[evidence_col].apply(get_xm_list)

    print(f"Writing output to: {out_path}")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Output file shape: {df.shape}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_barcodes_step1.py <bam_file> <fusion_tsv> <output_file>")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    tsv_path = sys.argv[2]
    out_path = sys.argv[3]
    
    print(f"Processing BAM: {bam_path}")
    print(f"Processing TSV: {tsv_path}")
    print(f"Output file: {out_path}")
    
    # Check file sizes
    try:
        bam_size = os.path.getsize(bam_path) / (1024**3)  # GB
        tsv_size = os.path.getsize(tsv_path) / (1024**2)  # MB
        print(f"BAM file size: {bam_size:.2f} GB")
        print(f"TSV file size: {tsv_size:.2f} MB")
    except:
        pass
    
    try:
        # Step 1: Get required read IDs from fusion TSV
        required_read_ids = get_required_read_ids(tsv_path)
        
        # Step 2: Extract tags only for required reads
        tag_dict = extract_tags_from_bam_selective(bam_path, required_read_ids)
        
        # Step 3: Update fusion TSV with tags
        update_fusion_tsv(tsv_path, out_path, tag_dict)
        
        print("Processing completed successfully.")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
