import pandas as pd
import os
import sys

def main():
    if len(sys.argv) != 3:
        print("Usage: python split_fusions_by_patient.py <input_tsv> <base_output_dir>")
        sys.exit(1)
    
    infile = sys.argv[1]
    base_outdir = sys.argv[2]

    # --- Extract pool_id from filename ---
    pool_id = os.path.basename(infile).split(".fusion_calls")[0]

    # make pool-specific directory
    outdir = os.path.join(base_outdir, pool_id)
    os.makedirs(outdir, exist_ok=True)

    # read TSV
    df = pd.read_csv(infile, sep="\t")

    # drop FFPM-related columns
    drop_cols = ["LR_FFPM", "max_LR_FFPMfrac_dom_iso", "above_frac_dom_iso"]
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])

    rows = []
    for _, row in df.iterrows():
        patients = str(row["assignment"]).split(",")
        n_patients = len(patients)

        # split num_LR evenly
        split_num = row["num_LR"] // n_patients
        remainder = row["num_LR"] % n_patients

        # split LR_accessions
        lr_reads = str(row["LR_accessions"]).split(",")
        chunk_size = len(lr_reads) // n_patients
        extra = len(lr_reads) % n_patients

        start = 0
        for i, pat in enumerate(patients):
            this_num = split_num + (1 if i < remainder else 0)
            this_chunk_size = chunk_size + (1 if i < extra else 0)
            this_reads = lr_reads[start:start+this_chunk_size]
            start += this_chunk_size

            new_row = row.copy()
            new_row["assignment"] = pat
            new_row["num_LR"] = this_num
            new_row["LR_accessions"] = ",".join(this_reads)
            new_row["pool_id"] = pool_id  # add pool_id column
            rows.append(new_row)

    # build patient-specific files inside pool directory
    out = pd.DataFrame(rows)
    for pat, subdf in out.groupby("assignment"):
        fname = os.path.join(outdir, f"{pat}.tsv")
        subdf.to_csv(fname, sep="\t", index=False)
        print(f"wrote {fname}")

if __name__ == "__main__":
    main()
