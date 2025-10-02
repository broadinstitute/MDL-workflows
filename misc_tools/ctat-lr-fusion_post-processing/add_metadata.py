import pandas as pd
import sys
from pathlib import Path

def revcomp(seq: str) -> str:
    """Function to reverse complement a DNA sequence"""
    if pd.isna(seq) or seq == "None" or seq == "":
        return seq
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]

def process_fusion_with_metadata(fusion_file, metadata_file, output_file, pool_name):
    # Load metadata
    print(f"Loading metadata from: {metadata_file}")
    meta = pd.read_csv(metadata_file, sep="\t")
    
    # Normalize barcodes
    meta = meta.rename(columns={"barcode": "CB"}) if "barcode" in meta.columns else meta
    meta["CB"] = meta["CB"].astype(str)
    
    # Also store reverse complement version
    meta["CB_rev"] = meta["CB"].apply(revcomp)
    
    # Extract pool_id from pool_name (e.g., "NeoMet_Pool_012_1" -> "NeoMet_Pool_012_1")
    pool_id = pool_name
    print(f"Using pool_id = {pool_id}")
    
    # Restrict metadata to this pool
    meta_pool = meta[meta["pool_id"] == pool_id]
    print(f"Found {len(meta_pool)} metadata records for pool {pool_id}")
    
    if meta_pool.empty:
        print(f"Warning: No metadata found for pool_id '{pool_id}'")
        print(f"Available pool_ids in metadata: {meta['pool_id'].unique()}")
        # Fallback to using all metadata if no pool-specific data found
        meta_pool = meta
        print("Using all metadata as fallback")
    
    # All metadata columns except CB/CB_rev and pool_id
    meta_fields = [c for c in meta.columns if c not in ("CB", "CB_rev", "pool_id")]
    print(f"Metadata fields to add: {meta_fields}")
    
    # Load fusion data
    print(f"Loading fusion data from: {fusion_file}")
    fusion = pd.read_csv(fusion_file, sep="\t")
    print(f"Fusion data shape: {fusion.shape}")
    
    # Initialize new columns
    for f in meta_fields:
        fusion[f] = ""
    
    # Iterate rows and map barcodes
    matched_rows = 0
    for i, row in fusion.iterrows():
        if pd.isna(row["CB_list"]) or row["CB_list"] == "":
            continue
            
        cbs = str(row["CB_list"]).split(",")
        cbs = [cb.strip() for cb in cbs if cb.strip() != "None" and cb.strip() != ""]
        
        if not cbs:
            continue
        
        # Try direct barcode match within the pool
        meta_sub = meta_pool[meta_pool["CB"].isin(cbs)]
        
        # If no match, try reverse complement
        if meta_sub.empty:
            rc_cbs = [revcomp(cb) for cb in cbs]
            meta_sub = meta_pool[meta_pool["CB"].isin(rc_cbs)]
        
        if not meta_sub.empty:
            matched_rows += 1
            for f in meta_fields:
                vals = meta_sub[f].dropna().astype(str).unique()
                vals = [v for v in vals if v != "nan" and v != ""]
                if len(vals) > 0:
                    fusion.at[i, f] = ",".join(vals)
    
    print(f"Successfully matched {matched_rows} fusion rows with metadata")
    
    # Write updated file
    fusion.to_csv(output_file, sep="\t", index=False)
    print(f"Output written to: {output_file}")
    print(f"Final output shape: {fusion.shape}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python add_metadata_step2.py <fusion_with_tags_tsv> <metadata_tsv> <output_file> <pool_name>")
        sys.exit(1)
    
    fusion_file = sys.argv[1]
    metadata_file = sys.argv[2]
    output_file = sys.argv[3]
    pool_name = sys.argv[4]
    
    print(f"Processing fusion file: {fusion_file}")
    print(f"Using metadata file: {metadata_file}")
    print(f"Pool name: {pool_name}")
    print(f"Output file: {output_file}")
    
    try:
        process_fusion_with_metadata(fusion_file, metadata_file, output_file, pool_name)
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
