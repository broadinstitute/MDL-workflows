version 1.0

workflow mark_duplicates_picard {
  input {
    String sample_id
    File bam

    # Optional: provide BAI if you already have one. If omitted, Picard can still run,
    # but it will also create an index for the output BAM when CREATE_INDEX=true.
    File? bai

    # Picard settings
    Boolean remove_duplicates = false
    Boolean assume_sorted = true
    Boolean create_index = true
    String validation_stringency = "SILENT"

    # Runtime
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb = 200

    # Use a Picard container; update to whatever you host/trust.
    String docker = "broadinstitute/picard:3.1.1"
  }

  call picard_markduplicates {
    input:
      sample_id = sample_id,
      bam = bam,
      bai = bai,
      remove_duplicates = remove_duplicates,
      assume_sorted = assume_sorted,
      create_index = create_index,
      validation_stringency = validation_stringency,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      docker = docker
  }

  output {
    File marked_bam = picard_markduplicates.marked_bam
    File metrics    = picard_markduplicates.metrics
    File marked_bai = picard_markduplicates.marked_bai
  }
}

task picard_markduplicates {
  input {
    String sample_id
    File bam
    File? bai

    Boolean remove_duplicates
    Boolean assume_sorted
    Boolean create_index
    String validation_stringency

    Int cpu
    Int memory_gb
    Int disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    # stable local names
    ln -s "~{bam}" in.bam
    if [[ "~{if defined(bai) then '1' else ''}" == "1" ]]; then
      ln -s "~{bai}" in.bam.bai
    fi

    picard \
      MarkDuplicates \
      I=in.bam \
      O="~{sample_id}.markdup.bam" \
      M="~{sample_id}.markdup.metrics.txt" \
      REMOVE_DUPLICATES="~{remove_duplicates}" \
      ASSUME_SORTED="~{assume_sorted}" \
      CREATE_INDEX="~{create_index}" \
      VALIDATION_STRINGENCY="~{validation_stringency}"
  >>>

  output {
    File marked_bam = "~{sample_id}.markdup.bam"
    File metrics    = "~{sample_id}.markdup.metrics.txt"
    File marked_bai = "~{sample_id}.markdup.bam.bai"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}
