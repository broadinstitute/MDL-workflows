version 1.0

task fastp_trim_pe {
  input {
    String sample_id

    File read1_fastq
    File read2_fastq

    # fastp parameters
    Int threads = 4
    Int min_length = 30
    Int qualified_quality_phred = 20
    Boolean detect_adapter_for_pe = true

    # Terra resources
    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb = 100

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/fastp-trimmer"
  }

  command <<<
    set -euo pipefail

    fastp \
      --in1 "~{read1_fastq}" \
      --in2 "~{read2_fastq}" \
      --out1 "~{sample_id}.R1.trim.fastq.gz" \
      --out2 "~{sample_id}.R2.trim.fastq.gz" \
      --thread ~{threads} \
      --qualified_quality_phred ~{qualified_quality_phred} \
      --length_required ~{min_length} \
      ~{if detect_adapter_for_pe then "--detect_adapter_for_pe" else ""} \
      --html "~{sample_id}.fastp.html" \
      --json "~{sample_id}.fastp.json" \
      2> "~{sample_id}.fastp.stderr.log"
  >>>

  output {
    File trimmed_read1 = "~{sample_id}.R1.trim.fastq.gz"
    File trimmed_read2 = "~{sample_id}.R2.trim.fastq.gz"
    File fastp_html = "~{sample_id}.fastp.html"
    File fastp_json = "~{sample_id}.fastp.json"
    File fastp_stderr_log = "~{sample_id}.fastp.stderr.log"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}

workflow fastp_trim_workflow {
  input {
    String sample_id
    File read1_fastq
    File read2_fastq

    Int threads = 4
    Int min_length = 30
    Int qualified_quality_phred = 20
    Boolean detect_adapter_for_pe = true

    Int cpu = 4
    Int memory_gb = 8
    Int disk_gb = 100
  }

  call fastp_trim_pe {
    input:
      sample_id = sample_id,
      read1_fastq = read1_fastq,
      read2_fastq = read2_fastq,
      threads = threads,
      min_length = min_length,
      qualified_quality_phred = qualified_quality_phred,
      detect_adapter_for_pe = detect_adapter_for_pe,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb
  }

  output {
    File trimmed_read1 = fastp_trim_pe.trimmed_read1
    File trimmed_read2 = fastp_trim_pe.trimmed_read2
    File fastp_html = fastp_trim_pe.fastp_html
    File fastp_json = fastp_trim_pe.fastp_json
    File fastp_stderr_log = fastp_trim_pe.fastp_stderr_log
  }
}
