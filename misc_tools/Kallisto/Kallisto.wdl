version 1.0

task KallistoQuantPairedEndSingle {
  input {
    String sample_id
    File kallisto_index
    File r1_fastq
    File r2_fastq

    Int threads = 8
    Int bootstraps = 0
    String strandedness = "unstranded"  # "unstranded" | "fr-stranded" | "rf-stranded"
    String extra_args = ""              # e.g. "--bias"
  }

  command <<<
    set -euo pipefail

    outdir="kallisto_out"
    mkdir -p "$outdir"

    STRAND_OPT=""
    if [[ "~{strandedness}" == "fr-stranded" ]]; then
      STRAND_OPT="--fr-stranded"
    elif [[ "~{strandedness}" == "rf-stranded" ]]; then
      STRAND_OPT="--rf-stranded"
    fi

    BS_OPT=""
    if [[ "~{bootstraps}" -gt 0 ]]; then
      BS_OPT="-b ~{bootstraps}"
    fi

    kallisto quant \
      -i "~{kallisto_index}" \
      -o "$outdir" \
      -t ~{threads} \
      $BS_OPT \
      $STRAND_OPT \
      ~{extra_args} \
      "~{r1_fastq}" "~{r2_fastq}" \
      2>&1 | tee "~{sample_id}.kallisto.log"

    # Rename outputs to include sample_id
    mv "$outdir/abundance.tsv" "$outdir/~{sample_id}.abundance.tsv"
    mv "$outdir/abundance.h5"  "$outdir/~{sample_id}.abundance.h5"
    mv "$outdir/run_info.json" "$outdir/~{sample_id}.run_info.json"

    tar -czf "~{sample_id}.kallisto_out.tgz" -C "$outdir" .
  >>>

  output {
    File abundance_tsv = "kallisto_out/~{sample_id}.abundance.tsv"
    File abundance_h5  = "kallisto_out/~{sample_id}.abundance.h5"
    File run_info_json = "kallisto_out/~{sample_id}.run_info.json"
    File log           = "~{sample_id}.kallisto.log"
    File output_tar    = "~{sample_id}.kallisto_out.tgz"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/methods-dev-lab/kallisto/kallisto"
    cpu: threads
    memory: "16G"
    disks: "local-disk 80 HDD"
  }
}

