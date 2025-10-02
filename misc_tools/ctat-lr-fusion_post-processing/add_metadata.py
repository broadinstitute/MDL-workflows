import pandas as pd
import sys

def revcomp(seq: str) -> str:
    """Function to reverse complement a DNA sequence"""
    if pd.isna(seq) or seq == "None" or seq == "":
        return seq
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]

def process_fusion_with_metadata(fusion_file, metadata_file, output_file):
    # Load metadata
    print(f"Loading metadata from: {metadata_file}")
    meta = pd.read_csv(metadata_file, sep="\t")
    
    # Normalize barcodes
    meta = meta.rename(columns={"barcode": "CB"}) if "barcode" in meta.columns else meta
    meta["CB"] = meta["CB"].astype(str)
    
    # Also store reverse complement version
    meta["CB_rev"] = meta["CB"].apply(revcomp)
    
    # All metadata columns except CB/CB_rev
    meta_fields = [c for c in meta.columns if c not in ("CB", "CB_rev")]
    
    # Load fusion data
    print(f"Loading fusion data from: {fusion_file}")
    fusion = pd.read_csv(fusion_file, sep="\t")
    
    # Initialize new columns
    for f in meta_fields:
        fusion[f] = ""
    
    # Iterate rows and map barcodes
    for i, row in fusion.iterrows():
        if pd.isna(row["CB_list"]) or row["CB_list"] == "":
            continue
            
        cbs = str(row["CB_list"]).split(",")
        cbs = [cb.strip() for cb in cbs if cb.strip() != "None" and cb.strip() != ""]
        
        if not cbs:
            continue
        
        # Check both direct and reverse complement matches
        meta_sub = meta[meta["CB"].isin(cbs)]
        if meta_sub.empty:
            rc_cbs = [revcomp(cb) for cb in cbs]
            meta_sub = meta[meta["CB"].isin(rc_cbs)]
        
        if not meta_sub.empty:
            for f in meta_fields:
                vals = meta_sub[f].dropna().astype(str).unique()
                vals = [v for v in vals if v != "nan" and v != ""]
                if len(vals) > 0:
                    fusion.at[i, f] = ",".join(vals)
    
    # Write updated file
    fusion.to_csv(output_file, sep="\t", index=False)
    print(f"Output written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python add_metadata_step2.py <fusion_with_tags_tsv> <metadata_tsv> <output_file>")
        sys.exit(1)
    
    fusion_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_file = sys.argv[3]
    
    process_fusion_with_metadata(fusion_file, metadata_file, output_file)
