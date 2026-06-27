version 1.0

workflow rseqc_gene_body_coverage {
  input {
    String analysis_id
    Array[String] sample_ids
    Array[File] bams
    Array[File] bais
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
      analysis_id = analysis_id,
      sample_ids = sample_ids,
      bams = bams,
      bais = bais,
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
    Array[File] heatmap_plots = run_gene_body_coverage.heatmap_plots
    File r_script = run_gene_body_coverage.r_script
    File sample_manifest = run_gene_body_coverage.sample_manifest
    File stdout_log = run_gene_body_coverage.stdout_log
    File stderr_log = run_gene_body_coverage.stderr_log
    File rseqc_log = run_gene_body_coverage.rseqc_log
    File output_tar_gz = run_gene_body_coverage.output_tar_gz
  }
}

task run_gene_body_coverage {
  input {
    String analysis_id
    Array[String] sample_ids
    Array[File] bams
    Array[File] bais
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

    MIN_MRNA_LENGTH="~{min_mrna_length}"
    OUTPUT_FORMAT="~{output_format}"
    REFGENE_BED="~{refgene_bed}"

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

    cat > sample_ids.list <<'EOF'
~{sep="\n" sample_ids}
EOF

    cat > bams.list <<'EOF'
~{sep="\n" bams}
EOF

    cat > bais.list <<'EOF'
~{sep="\n" bais}
EOF

    sample_count=$(wc -l < sample_ids.list)
    bam_count=$(wc -l < bams.list)
    bai_count=$(wc -l < bais.list)

    if [[ "${sample_count}" -eq 0 ]]; then
      echo "ERROR: At least one sample is required." >&2
      exit 1
    fi

    lengths_match=1
    if (( sample_count != bam_count )); then
      lengths_match=0
    fi
    if (( sample_count != bai_count )); then
      lengths_match=0
    fi

    if (( lengths_match == 0 )); then
      echo "ERROR: sample_ids, bams, and bais must have the same length." >&2
      echo "sample_ids=${sample_count} bams=${bam_count} bais=${bai_count}" >&2
      exit 1
    fi

    if [[ "${REFGENE_BED}" == *.gz ]]; then
      gzip -dc "${REFGENE_BED}" > refgene.bed
    else
      ln -s "${REFGENE_BED}" refgene.bed
    fi

    mkdir -p staged_bams
    : > bam_paths.list
    printf "sample_id\tstaged_bam\n" > sample_manifest.tsv

    declare -A seen_names=()
    index=0
    while IFS= read -r sample_id && IFS= read -r bam <&3 && IFS= read -r bai <&4; do
      index=$((index + 1))
      safe_name=$(printf "%s" "${sample_id}" | sed -E 's/[^A-Za-z0-9_.]+/_/g; s/^_+//; s/_+$//')
      if [[ -z "${safe_name}" ]]; then
        safe_name="sample_${index}"
      fi

      if [[ -n "${seen_names[${safe_name}]:-}" ]]; then
        safe_name="${safe_name}_${index}"
      fi
      seen_names["${safe_name}"]=1

      ln -s "${bam}" "staged_bams/${safe_name}.bam"
      ln -s "${bai}" "staged_bams/${safe_name}.bam.bai"
      printf "staged_bams/%s.bam\n" "${safe_name}" >> bam_paths.list
      printf "%s\t%s.bam\n" "${sample_id}" "${safe_name}" >> sample_manifest.tsv
    done < sample_ids.list 3< bams.list 4< bais.list

    geneBody_coverage.py \
      --input bam_paths.list \
      --refgene refgene.bed \
      --minimum_length "${MIN_MRNA_LENGTH}" \
      --format "${OUTPUT_FORMAT}" \
      --out-prefix "~{analysis_id}" \
      > "~{analysis_id}.geneBodyCoverage.stdout.log" \
      2> "~{analysis_id}.geneBodyCoverage.stderr.log"

    test -s "~{analysis_id}.geneBodyCoverage.txt"
    test -s "~{analysis_id}.geneBodyCoverage.curves.~{output_format}"
    test -s "~{analysis_id}.geneBodyCoverage.r"
    if [[ "${sample_count}" -ge 3 ]]; then
      test -s "~{analysis_id}.geneBodyCoverage.heatMap.~{output_format}"
    fi

    tar_files=(
      "~{analysis_id}.geneBodyCoverage.txt" \
      "~{analysis_id}.geneBodyCoverage.r" \
      "~{analysis_id}.geneBodyCoverage.curves.~{output_format}" \
      "~{analysis_id}.geneBodyCoverage.stdout.log" \
      "~{analysis_id}.geneBodyCoverage.stderr.log" \
      sample_manifest.tsv \
      log.txt
    )
    if [[ "${sample_count}" -ge 3 ]]; then
      tar_files+=("~{analysis_id}.geneBodyCoverage.heatMap.~{output_format}")
    fi

    tar -czf "~{analysis_id}.geneBodyCoverage.outputs.tar.gz" "${tar_files[@]}"
  >>>

  output {
    File coverage_table = "~{analysis_id}.geneBodyCoverage.txt"
    File curves_plot = "~{analysis_id}.geneBodyCoverage.curves.~{output_format}"
    Array[File] heatmap_plots = glob("~{analysis_id}.geneBodyCoverage.heatMap.~{output_format}")
    File r_script = "~{analysis_id}.geneBodyCoverage.r"
    File sample_manifest = "sample_manifest.tsv"
    File stdout_log = "~{analysis_id}.geneBodyCoverage.stdout.log"
    File stderr_log = "~{analysis_id}.geneBodyCoverage.stderr.log"
    File rseqc_log = "log.txt"
    File output_tar_gz = "~{analysis_id}.geneBodyCoverage.outputs.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible
  }
}
