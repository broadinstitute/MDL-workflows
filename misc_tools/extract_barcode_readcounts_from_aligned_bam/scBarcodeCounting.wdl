version 1.0

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
      --processes ~{cpu} \
      --threads 2 \
      ~{skip_umi_arg} \
      -o ~{sample_id}.counts.tsv
  >>>

  output {
    File counts = "~{sample_id}.counts.tsv"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "32G"
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
