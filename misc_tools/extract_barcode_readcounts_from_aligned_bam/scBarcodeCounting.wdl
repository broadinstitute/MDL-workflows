version 1.0

task CountUMIs {
  input {
    File bam
    File bai
    File process_barcodes_py
    String sample_id
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  command <<<
    set -euox pipefail

    cp ~{process_barcodes_py} process_barcodes_py.py

    python3 process_barcodes_py.py \
      ~{bam} \
      --sample_id ~{sample_id} \
      --cb_tag ~{cb_tag} \
      --umi_tag ~{umi_tag} \
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

task MergeCounts {
  input {
    Array[File] count_files
  }

  command <<<
    set -euo pipefail

    # Just concatenate all per-sample TSVs into one
    head -n 1 ~{count_files[0]} > merged_counts.tsv
    for f in ~{sep=" " count_files}; do
      tail -n +2 "$f" >> merged_counts.tsv
    done
  >>>

  output {
    File merged = "merged_counts.tsv"
  }

  runtime {
    docker: "ubuntu:20.04"
    cpu: 1
    memory: "2G"
    disks: "local-disk 10 SSD"
  }
}

workflow CountUMIsWorkflow {
  input {
    Array[File] bams
    Array[File] bais
    Array[String] sample_ids
    File process_barcodes_py
    String cb_tag = "CB"
    String umi_tag = "XM"
    Int threads = 8
  }

  scatter (i in range(length(bams))) {
    call CountUMIs {
      input:
        bam = bams[i],
        bai = bais[i],
        process_barcodes_py = process_barcodes_py,
        sample_id = sample_ids[i],
        cb_tag = cb_tag,
        umi_tag = umi_tag,
        threads = threads
    }
  }

  call MergeCounts {
    input:
      count_files = CountUMIs.counts
  }

  output {
    Array[File] all_counts = CountUMIs.counts
    File merged_counts = MergeCounts.merged
  }
}
