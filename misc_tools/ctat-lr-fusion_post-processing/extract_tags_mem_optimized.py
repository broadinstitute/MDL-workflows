import pysam
import pandas as pd
from pathlib import Path
import sys
from datetime import datetime

def setup_logging(log_path):
    """
    Setup logging to file.
    """
    class Logger:
        def __init__(self, filepath):
            self.terminal = sys.stdout
            self.log = open(filepath, 'w')
        
        def write(self, message):
            self.terminal.write(message)
            self.log.write(message)
            self.log.flush()  # Ensure immediate write
        
        def flush(self):
            self.terminal.flush()
            self.log.flush()
    
    return Logger(log_path)

def log_message(msg, log_file=None):
    """
    Print message and write to log file if provided.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    message = f"[{timestamp}] {msg}"
    print(message, flush=True)
    if log_file:
        log_file.write(message + "\n")
        log_file.flush()

def get_required_read_ids(tsv_path, log_file=None):
    """
    Extract all unique read IDs from the fusion TSV evidence column.
    Returns a set of read IDs that we need to look up in the BAM.
    """
    log_message(f"Reading TSV file: {tsv_path}", log_file)
    df = pd.read_csv(tsv_path, sep="\t")
    evidence_col = df.columns[9]
    log_message(f"Evidence column: {evidence_col}", log_file)
    
    read_ids = set()
    for reads_str in df[evidence_col]:
        if pd.notna(reads_str):
            reads = str(reads_str).split(",")
            read_ids.update(reads)
    
    log_message(f"Found {len(read_ids)} unique read IDs to look up", log_file)
    return read_ids

def extract_tags_from_bam(bam_path, required_read_ids, log_file=None):
    """
    Returns a dict mapping read_id -> (XM, CB).
    Only indexes reads that are in required_read_ids.
    """
    log_message(f"Processing BAM file: {bam_path}", log_file)
    log_message(f"Looking for {len(required_read_ids)} reads", log_file)
    
    # Initialize dict with all required read IDs
    tag_dict = {read_id: (None, None) for read_id in required_read_ids}
    
    found_count = 0
    total_reads_scanned = 0
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            total_reads_scanned += 1
            
            if total_reads_scanned % 1000000 == 0:
                log_message(f"Scanned {total_reads_scanned:,} reads, found {found_count}/{len(required_read_ids)} targets", log_file)
            
            read_id = read.query_name
            # Only process if this read is in our required set
            if read_id in tag_dict:
                xm = read.get_tag("XM") if read.has_tag("XM") else None
                cb = read.get_tag("CB") if read.has_tag("CB") else None
                tag_dict[read_id] = (xm, cb)
                found_count += 1
                
                # Early exit if we've found all required reads
                if found_count == len(required_read_ids):
                    log_message(f"Found all {found_count} required reads after scanning {total_reads_scanned:,} total reads", log_file)
                    break
    
    log_message(f"BAM processing complete. Found {found_count}/{len(required_read_ids)} reads", log_file)
    
    if found_count < len(required_read_ids):
        missing = len(required_read_ids) - found_count
        log_message(f"WARNING: {missing} read IDs from TSV were not found in BAM", log_file)
    
    return tag_dict

def update_fusion_tsv(tsv_path, out_path, tag_dict, log_file=None):
    log_message(f"Updating fusion TSV", log_file)
    df = pd.read_csv(tsv_path, sep="\t")
    evidence_col = df.columns[9]

    def get_cb_list(reads_str):
        if pd.isna(reads_str):
            return ""
        reads = str(reads_str).split(",")
        return ",".join([str(tag_dict.get(r, (None, None))[1]) if tag_dict.get(r, (None, None))[1] is not None else "None" for r in reads])

    def get_xm_list(reads_str):
        if pd.isna(reads_str):
            return ""
        reads = str(reads_str).split(",")
        return ",".join([str(tag_dict.get(r, (None, None))[0]) if tag_dict.get(r, (None, None))[0] is not None else "None" for r in reads])

    df["CB_list"] = df[evidence_col].apply(get_cb_list)
    df["XM_list"] = df[evidence_col].apply(get_xm_list)

    log_message(f"Writing output to: {out_path}", log_file)
    df.to_csv(out_path, sep="\t", index=False)
    log_message(f"Output written successfully", log_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_barcodes_step1.py <bam_file> <fusion_tsv> <output_file>")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    tsv_path = sys.argv[2]
    out_path = sys.argv[3]
    
    # Create log file in same directory as output
    log_path = str(Path(out_path).parent / "extract_barcodes.log")
    
    # Open log file
    with open(log_path, 'w') as log_file:
        log_message("="*60, log_file)
        log_message("Starting barcode extraction process", log_file)
        log_message("="*60, log_file)
        log_message(f"BAM file: {bam_path}", log_file)
        log_message(f"TSV file: {tsv_path}", log_file)
        log_message(f"Output file: {out_path}", log_file)
        log_message(f"Log file: {log_path}", log_file)
        log_message("="*60, log_file)
        
        try:
            required_read_ids = get_required_read_ids(tsv_path, log_file)
            tag_dict = extract_tags_from_bam(bam_path, required_read_ids, log_file)
            update_fusion_tsv(tsv_path, out_path, tag_dict, log_file)
            
            log_message("="*60, log_file)
            log_message("Processing complete!", log_file)
            log_message("="*60, log_file)
        except Exception as e:
            log_message("="*60, log_file)
            log_message(f"ERROR: {str(e)}", log_file)
            log_message("="*60, log_file)
            import traceback
            log_file.write(traceback.format_exc())
            raise
    
    print(f"\n📄 Log file written to: {log_path}")
