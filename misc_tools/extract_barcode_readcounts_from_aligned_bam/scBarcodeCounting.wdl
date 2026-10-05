version 1.0

task CountUMI {
  input {
    File bam
    File bai
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.0"
  }

  command <<<
    set -euo pipefail

    process_barcodes.py \
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
    docker: docker_image
    cpu: 4
    memory: "32G"
    disks: "local-disk 200 SSD"
  }
}

workflow CountOneBAMWorkflow {
  input {
    File bam
    File bai
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.0"
  }

  call CountUMI {
    input:
      bam = bam,
      bai = bai,
      sample_id = sample_id,
      cb_tag = cb_tag,
      umi_tag = umi_tag,
      docker_image = docker_image
  }

  output {
    File counts = CountUMI.counts
  }
}
