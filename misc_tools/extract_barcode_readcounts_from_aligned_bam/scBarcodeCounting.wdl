version 1.0

task CountUMI {
  input {
    File bam
    File bai
    File process_barcodes_py # Python script
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
  }

  command <<<
    set -euox pipefail

    cp ~{process_barcodes_py} process_barcodes_py.py

    python3 process_barcodes_py.py \
      ~{bam} \
      --sample_id ~{sample_id} \
      --cb_tag ~{cb_tag} \
      --umi_tag ~{umi_tag} \
      -o ~{sample_id}.counts.tsv
  >>>

  output {
    File counts = "~{sample_id}.counts.tsv"
  }

  runtime {
    docker: "python:3.10-slim"
    cpu: 4
    memory: "16G"
    disks: "local-disk 100 SSD"
  }
}

workflow CountOneBAMWorkflow {
  input {
    File bam
    File bai
    File process_barcodes_py
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
  }

  call CountUMI {
    input:
      bam = bam,
      bai = bai,
      process_barcodes_py = process_barcodes_py,
      sample_id = sample_id,
      cb_tag = cb_tag,
      umi_tag = umi_tag
  }

  output {
    File counts = CountUMI.counts
  }
}
