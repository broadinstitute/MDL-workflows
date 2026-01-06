version 1.0

workflow rnaseqc2 {
  input {
    String sample_id
    File bam
    File? bai

    # IMPORTANT: RNA-SeQC 2 expects a "collapsed" GTF (see README).
    File collapsed_gtf

    # Optional: if you want fragment size metrics, provide a BED of non-overlapping exons
    File? non_overlapping_exons_bed

    # Optional: only needed for CRAM (or if you want reference-based mismatch filtering)
    File? reference_fasta

    # Optional strandedness mode: "RF", "FR", "rf", "fr"
    String? stranded

    # Optional knobs
    Int mapping_quality = 255
    Int detection_threshold = 5
    Int threads = 4
    Int memory_gb = 16

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/rnaseqc2/rnaseqc2"
  }

  call ensure_bam_index {
    input:
      bam = bam,
      bai = bai,
      threads = threads,
      memory_gb = memory_gb,
      docker = docker
  }

  call run_rnaseqc {
    input:
      sample_id = sample_id,
      bam = ensure_bam_index.bam_out,
      bai = ensure_bam_index.bai_out,
      collapsed_gtf = collapsed_gtf,
      non_overlapping_exons_bed = non_overlapping_exons_bed,
      reference_fasta = reference_fasta,
      stranded = stranded,
      mapping_quality = mapping_quality,
      detection_threshold = detection_threshold,
      threads = threads,
      memory_gb = memory_gb,
      docker = docker
  }

  output {
    File metrics_tsv       = run_rnaseqc.metrics_tsv
    File gene_reads_gct    = run_rnaseqc.gene_reads_gct
    File gene_tpm_gct      = run_rnaseqc.gene_tpm_gct
    File exon_reads_gct    = run_rnaseqc.exon_reads_gct
    File output_dir_tar_gz = run_rnaseqc.output_dir_tar_gz
  }
}

task ensure_bam_index {
  input {
    File bam
    File? bai
    Int threads
    Int memory_gb
    String docker
  }

  command <<<
    set -euo pipefail

    ln -s "~{bam}" sample.bam

    if [[ "~{if defined(bai) then '1' else ''}" == "1" ]]; then
      ln -s "~{bai}" sample.bam.bai
    else
      # RNA-SeQC reads BAMs; it doesn't require an index for all operations,
      # but providing one is generally safer and helps if downstream tooling expects it.
      # Use samtools if present; otherwise, fail with a clear message.
      if command -v samtools >/dev/null 2>&1; then
        samtools index -@ "~{threads}" sample.bam
      else
        echo "ERROR: bai not provided and samtools not available in container to build one." >&2
        exit 1
      fi
    fi
  >>>

  output {
    File bam_out = "sample.bam"
    File bai_out = "sample.bam.bai"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GiB"
  }
}

task run_rnaseqc {
  input {
    String sample_id
    File bam
    File bai
    File collapsed_gtf
    File? non_overlapping_exons_bed
    File? reference_fasta
    String? stranded

    Int mapping_quality
    Int detection_threshold

    Int threads
    Int memory_gb
    String docker
  }

  command <<<
    set -euo pipefail

    mkdir -p out

    # Keep filenames stable
    ln -s "~{bam}" "~{sample_id}.bam"
    ln -s "~{bai}" "~{sample_id}.bam.bai"

    # Build optional flags
    EXTRA_ARGS=()
    EXTRA_ARGS+=(--mapping-quality "~{mapping_quality}")
    EXTRA_ARGS+=(--detection-threshold "~{detection_threshold}")

    if [[ -n "~{if defined(non_overlapping_exons_bed) then '1' else ''}" ]]; then
      ln -s "~{non_overlapping_exons_bed}" exons.bed
      EXTRA_ARGS+=(--bed exons.bed)
    fi

    if [[ -n "~{if defined(reference_fasta) then '1' else ''}" ]]; then
      ln -s "~{reference_fasta}" ref.fa
      EXTRA_ARGS+=(--fasta ref.fa)
    fi

    if [[ -n "~{if defined(stranded) then stranded else ''}" ]]; then
      EXTRA_ARGS+=(--stranded "~{stranded}")
    fi

    # Run RNA-SeQC 2
    # CLI: rnaseqc [OPTIONS] gtf bam output_dir
    rnaseqc "${EXTRA_ARGS[@]}" \
      "~{collapsed_gtf}" \
      "~{sample_id}.bam" \
      out \
      --sample "~{sample_id}" \
      -v

    # Bundle full output directory for convenience in Terra
    tar -czf "~{sample_id}.rnaseqc2.outputs.tar.gz" -C out .
  >>>

  output {
    File metrics_tsv    = "out/~{sample_id}.metrics.tsv"
    File exon_reads_gct = "out/~{sample_id}.exon_reads.gct"
    File gene_reads_gct = "out/~{sample_id}.gene_reads.gct"
    File gene_tpm_gct   = "out/~{sample_id}.gene_tpm.gct"

    File output_dir_tar_gz = "~{sample_id}.rnaseqc2.outputs.tar.gz"
  }

  runtime {
    docker: docker
    cpu: threads
    memory: "~{memory_gb} GiB"
  }
}
