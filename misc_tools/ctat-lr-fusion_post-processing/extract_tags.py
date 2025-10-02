import pysam
import pandas as pd
from pathlib import Path
import sys

def extract_tags_from_bam(bam_path):
    """
    Returns a dict mapping read_id -> (XM, CB).
    """
    tag_dict = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            read_id = read.query_name
            xm = read.get_tag("XM") if read.has_tag("XM") else None
            cb = read.get_tag("CB") if read.has_tag("CB") else None
            tag_dict[read_id] = (xm, cb)
    return tag_dict

def update_fusion_tsv(tsv_path, out_path, tag_dict):
    df = pd.read_csv(tsv_path, sep="\t")

    # assuming evidence reads are in the 10th column (index 9)
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

    df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_barcodes_step1.py <bam_file> <fusion_tsv> <output_file>")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    tsv_path = sys.argv[2]
    out_path = sys.argv[3]
    
    print(f"Processing BAM: {bam_path}")
    print(f"Processing TSV: {tsv_path}")
    
    tag_dict = extract_tags_from_bam(bam_path)
    update_fusion_tsv(tsv_path, out_path, tag_dict)
    
    print("✅ Processing complete.")
