version 1.0

task CountUMIs {
  input {
    File bam
    File bai
    File counting_script   # <-- Python script provided as input
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  command <<<
    set -euox pipefail

    cp ~{counting_script} process_barcodes.py

    python3 ./process_barcodes.py \
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
    File counting_script
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  scatter (i in range(length(bams))) {
    call CountUMIs {
      input:
        bam = bams[i],
        bai = bais[i],
        counting_script = counting_script,
        cb_tag = cb_tag,
        umi_tag = umi_tag,
        threads = threads
    }
  }

  output {
    Array[File] all_counts = CountUMIs.counts
  }
}
