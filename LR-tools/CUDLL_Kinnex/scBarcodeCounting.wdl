version 1.1

task CountUMI {
  input {
    File bam
    File bai
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
    Boolean skip_umi = false
    Int cpu = 8
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.1"
  }

  Int diskGB = ceil(size(bam, "GB") * 1.2 + size(bai, "GB")) + 20

  # Memory scales with input size: peak was 2.03 GiB for a ~20 GB BAM (one data point), so
  # 0.3 GB per GB of BAM + 4 GB gives ~5x headroom there. Floor of 8 GB.
  Int scaled_memory_gb = ceil(size(bam, "GB") * 0.3) + 4
  Int task_memory_gb = if scaled_memory_gb > 8 then scaled_memory_gb else 8

  # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
  Int cpu_tier = if cpu <= 4 then 4
    else if cpu <= 8 then 8
    else if cpu <= 16 then 16
    else if cpu <= 30 then 30
    else if cpu <= 60 then 60
    else if cpu <= 90 then 90
    else if cpu <= 180 then 180
    else 360
  Int mem_tier = if task_memory_gb <= 32 then 4
    else if task_memory_gb <= 64 then 8
    else if task_memory_gb <= 128 then 16
    else if task_memory_gb <= 240 then 30
    else if task_memory_gb <= 480 then 60
    else if task_memory_gb <= 720 then 90
    else if task_memory_gb <= 1440 then 180
    else 360
  Int effective_cpu = if cpu_tier >= mem_tier then cpu_tier else mem_tier
  # c3d-highcpu RAM is non-uniform per tier; use exact values.
  Int highcpu_ram = if effective_cpu == 4 then 8
    else if effective_cpu == 8 then 16
    else if effective_cpu == 16 then 32
    else if effective_cpu == 30 then 59
    else if effective_cpu == 60 then 118
    else if effective_cpu == 90 then 177
    else if effective_cpu == 180 then 354
    else 708
  String machine_type = if task_memory_gb <= highcpu_ram
    then "c3d-highcpu-${effective_cpu}"
    else if task_memory_gb <= effective_cpu * 4
    then "c3d-standard-${effective_cpu}"
    else "c3d-highmem-${effective_cpu}"

  String skip_umi_arg = if skip_umi then "--skip_umi" else ""

  command <<<
    set -euo pipefail

    # pysam needs the index next to the BAM as <bam>.bai
    ln -s ~{bam} input.bam
    ln -s ~{bai} input.bam.bai

    # One worker per CPU, each taking a share of the contigs (balanced from the .bai).
    process_barcodes.py input.bam \
      --sample_id ~{sample_id} \
      --cb_tag ~{cb_tag} \
      --umi_tag ~{umi_tag} \
      --processes ~{effective_cpu} \
      --threads 2 \
      ~{skip_umi_arg} \
      -o ~{sample_id}.counts.tsv
  >>>

  output {
    File counts = "~{sample_id}.counts.tsv"
  }

  # GCP Batch rejects predefinedMachineType + explicit cpu/memory together; omit them.
  runtime {
    docker: docker_image
    predefinedMachineType: "~{machine_type}"
    disks: "local-disk ~{diskGB} SSD"
    preemptible: 2
  }
}

# Sums reads per barcode across all BAMs into a barcode/post_count table, the default
# input of Shard_Bams_By_Barcode_Group (counts_column=post_count, barcode_column=barcode).
task MergeCounts {
  input {
    Array[File] counts
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.1"
  }

  Int diskGB = ceil(size(counts, "GB") * 3) + 10

  command <<<
    set -euo pipefail

    merge_barcode_counts.py ~{sep=' ' counts} -o merged_counts.tsv
  >>>

  output {
    File merged_counts = "merged_counts.tsv"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "8G"
    disks: "local-disk ~{diskGB} HDD"
    preemptible: 2
  }
}

workflow CountBarcodesWorkflow {
  input {
    Array[File] bams
    Array[File] bais
    String cb_tag = "CB"
    String umi_tag = "XM"
    Boolean skip_umi = false
    Int cpu = 8
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.1"
  }

  scatter (i in range(length(bams))) {
    # sample_id = BAM file name without the .bam extension
    call CountUMI {
      input:
        bam = bams[i],
        bai = bais[i],
        sample_id = sub(basename(bams[i]), "\\.bam$", ""),
        cb_tag = cb_tag,
        umi_tag = umi_tag,
        skip_umi = skip_umi,
        cpu = cpu,
        docker_image = docker_image
    }
  }

  call MergeCounts {
    input:
      counts = CountUMI.counts,
      docker_image = docker_image
  }

  output {
    Array[File] counts = CountUMI.counts          # per BAM: sample_id, cell_barcode, reads, umis
    File merged_counts = MergeCounts.merged_counts # barcode, post_count (all BAMs summed)
  }
}

