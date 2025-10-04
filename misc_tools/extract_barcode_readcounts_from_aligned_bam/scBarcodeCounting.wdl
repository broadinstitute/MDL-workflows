version 1.0

task CountUMIs {
  input {
    File bam
    File bai
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  command <<<
    set -euo pipefail

    python3 count_barcodes.py \
      ~{bam} \
      --cb_tag ~{cb_tag} \
      --umi_tag ~{umi_tag} \
      --threads ~{threads} \
      -o counts.tsv
  >>>

  output {
    File counts = "counts.tsv"
  }

  runtime {
    docker: "python:3.10-slim"
    cpu: threads
    memory: "16G"
    disks: "local-disk 100 SSD"
  }
}

workflow CountUMIsWorkflow {
  input {
    Array[File] bams
    Array[File] bais
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  scatter (i in range(length(bams))) {
    call CountUMIs {
      input:
        bam = bams[i],
        bai = bais[i],
        cb_tag = cb_tag,
        umi_tag = umi_tag,
        threads = threads
    }
  }

  output {
    Array[File] all_counts = CountUMIs.counts
  }
}
