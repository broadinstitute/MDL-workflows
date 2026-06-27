version 1.0

workflow rseqc_gene_body_coverage {
  input {
    String sample_id
    File bam
    File bai
    File refgene_bed

    Int min_mrna_length = 100
    String output_format = "pdf"

    Int cpu = 2
    Int memory_gb = 16
    Int disk_gb = 100
    Int preemptible = 1

    String docker = "quay.io/biocontainers/rseqc:5.0.4--pyhdfd78af_1"
  }

  call run_gene_body_coverage {
    input:
      sample_id = sample_id,
      bam = bam,
      bai = bai,
      refgene_bed = refgene_bed,
      min_mrna_length = min_mrna_length,
      output_format = output_format,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      preemptible = preemptible,
      docker = docker
  }

  output {
    File coverage_table = run_gene_body_coverage.coverage_table
    File curves_plot = run_gene_body_coverage.curves_plot
    File r_script = run_gene_body_coverage.r_script
    File stdout_log = run_gene_body_coverage.stdout_log
    File stderr_log = run_gene_body_coverage.stderr_log
    File rseqc_log = run_gene_body_coverage.rseqc_log
    File multiqc_files_tar_gz = run_gene_body_coverage.multiqc_files_tar_gz
    File output_tar_gz = run_gene_body_coverage.output_tar_gz
  }
}

task run_gene_body_coverage {
  input {
    String sample_id
    File bam
    File bai
    File refgene_bed

    Int min_mrna_length
    String output_format

    Int cpu
    Int memory_gb
    Int disk_gb
    Int preemptible

    String docker
  }

  command <<<
    set -euo pipefail

    SAMPLE_ID="~{sample_id}"
    MIN_MRNA_LENGTH="~{min_mrna_length}"
    OUTPUT_FORMAT="~{output_format}"
    REFGENE_BED="~{refgene_bed}"

    if [[ -z "${SAMPLE_ID}" ]]; then
      echo "ERROR: sample_id cannot be empty." >&2
      exit 1
    fi

    if [[ ! "${SAMPLE_ID}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
      echo "ERROR: sample_id can only contain letters, numbers, '.', '_', and '-'." >&2
      exit 1
    fi

    if [[ "${MIN_MRNA_LENGTH}" -lt 100 ]]; then
      echo "ERROR: min_mrna_length must be >= 100 for geneBody_coverage.py." >&2
      exit 1
    fi

    case "${OUTPUT_FORMAT}" in
      pdf|png|jpeg) ;;
      *)
        echo "ERROR: output_format must be one of: pdf, png, jpeg." >&2
        exit 1
        ;;
    esac

    if [[ "${REFGENE_BED}" == *.gz ]]; then
      gzip -dc "${REFGENE_BED}" > refgene.bed
    else
      ln -s "${REFGENE_BED}" refgene.bed
    fi

    ln -s "~{bam}" "${SAMPLE_ID}.bam"
    ln -s "~{bai}" "${SAMPLE_ID}.bam.bai"

    geneBody_coverage.py \
      --input "${SAMPLE_ID}.bam" \
      --refgene refgene.bed \
      --minimum_length "${MIN_MRNA_LENGTH}" \
      --format "${OUTPUT_FORMAT}" \
      --out-prefix "${SAMPLE_ID}" \
      > "${SAMPLE_ID}.geneBodyCoverage.stdout.log" \
      2> "${SAMPLE_ID}.geneBodyCoverage.stderr.log"

    mv log.txt "${SAMPLE_ID}.geneBodyCoverage.rseqc.log"

    test -s "${SAMPLE_ID}.geneBodyCoverage.txt"
    test -s "${SAMPLE_ID}.geneBodyCoverage.curves.${OUTPUT_FORMAT}"
    test -s "${SAMPLE_ID}.geneBodyCoverage.r"
    test -s "${SAMPLE_ID}.geneBodyCoverage.rseqc.log"

    tar_files=(
      "${SAMPLE_ID}.geneBodyCoverage.txt"
      "${SAMPLE_ID}.geneBodyCoverage.r"
      "${SAMPLE_ID}.geneBodyCoverage.curves.${OUTPUT_FORMAT}"
      "${SAMPLE_ID}.geneBodyCoverage.stdout.log"
      "${SAMPLE_ID}.geneBodyCoverage.stderr.log"
      "${SAMPLE_ID}.geneBodyCoverage.rseqc.log"
    )

    tar -czf "${SAMPLE_ID}.geneBodyCoverage.outputs.tar.gz" "${tar_files[@]}"
    tar -czf "${SAMPLE_ID}.geneBodyCoverage.multiqc_files.tar.gz" \
      "${SAMPLE_ID}.geneBodyCoverage.txt"
  >>>

  output {
    File coverage_table = "~{sample_id}.geneBodyCoverage.txt"
    File curves_plot = "~{sample_id}.geneBodyCoverage.curves.~{output_format}"
    File r_script = "~{sample_id}.geneBodyCoverage.r"
    File stdout_log = "~{sample_id}.geneBodyCoverage.stdout.log"
    File stderr_log = "~{sample_id}.geneBodyCoverage.stderr.log"
    File rseqc_log = "~{sample_id}.geneBodyCoverage.rseqc.log"
    File multiqc_files_tar_gz = "~{sample_id}.geneBodyCoverage.multiqc_files.tar.gz"
    File output_tar_gz = "~{sample_id}.geneBodyCoverage.outputs.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }
}
