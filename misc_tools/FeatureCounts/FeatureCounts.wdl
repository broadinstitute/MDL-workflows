version 1.0

workflow featurecounts_workflow {
  input {
    String sample_id
    File bam
    File? bai
    File gtf

    # featureCounts options
    Boolean is_paired_end = true
    Int strandness = 0              # 0=unstranded, 1=stranded, 2=reverse-stranded
    String feature_type = "exon"    # -t
    String gene_attr = "gene_id"    # -g

    Int threads = 4
    String extra_args = ""          # any additional featureCounts flags

    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/featurecounts/featurecounts"
  }

  call run_featurecounts {
    input:
      sample_id = sample_id,
      bam = bam,
      bai = bai,
      gtf = gtf,
      is_paired_end = is_paired_end,
      strandness = strandness,
      feature_type = feature_type,
      gene_attr = gene_attr,
      threads = threads,
      extra_args = extra_args,
      cpu = cpu,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      docker = docker
  }

  output {
    File counts = run_featurecounts.counts
    File summary = run_featurecounts.summary
    File log = run_featurecounts.log
  }
}

task run_featurecounts {
  input {
    String sample_id
    File bam
    File? bai
    File gtf

    Boolean is_paired_end = true
    Int strandness = 0
    String feature_type = "exon"
    String gene_attr = "gene_id"

    Int threads = 4
    String extra_args = ""

    Int cpu = 4
    Int memory_gb = 32
    Int disk_gb = 100

    String docker
  }

  command <<<
    set -euo pipefail

    # Use a sample-named BAM in the working dir (nicer output names / multiqc)
    ln -s "~{bam}" "~{sample_id}.bam"
    if [[ -n "~{select_first([bai, ''])}" ]]; then
      ln -s "~{select_first([bai, ''])}" "~{sample_id}.bam.bai" || true
    fi

    # Localize GTF; if gzipped, decompress to a plain .gtf for safety/compatibility
    GTF_LOCAL="annotation.gtf"
    if [[ "~{gtf}" == *.gz ]]; then
      zcat "~{gtf}" > "${GTF_LOCAL}"
    else
      ln -s "~{gtf}" "${GTF_LOCAL}"
    fi

    OUT_PREFIX="~{sample_id}.featureCounts"
    OUT_TXT="${OUT_PREFIX}.txt"

    PAIRED_FLAG=""
    if [[ "~{is_paired_end}" == "true" ]]; then
      PAIRED_FLAG="-p"
    fi

    # Run featureCounts
    # Notes:
    #  -T threads
    #  -s strandness (0/1/2)
    #  -t feature_type (exon for gene-level counting)
    #  -g gene_attr (gene_id for Gencode)
    featureCounts \
      -T "~{threads}" \
      -a "${GTF_LOCAL}" \
      -o "${OUT_TXT}" \
      ${PAIRED_FLAG} \
      -s "~{strandness}" \
      -t "~{feature_type}" \
      -g "~{gene_attr}" \
      ~{extra_args} \
      "~{sample_id}.bam" \
      2> "~{sample_id}.featureCounts.stderr.log"

    # featureCounts writes summary as <out>.summary (i.e. OUT_TXT.summary)
    ls -lh
  >>>

  output {
    File counts = "~{sample_id}.featureCounts.txt"
    File summary = "~{sample_id}.featureCounts.txt.summary"
    File log = "~{sample_id}.featureCounts.stderr.log"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
  }
}
