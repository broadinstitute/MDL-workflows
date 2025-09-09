version 1.0


workflow BamErrorRates_wf {
  input {
    File bam_file
    String sample_id
  }

  call BamErrorRates {
    input:
      input_bam = bam_file,
      sample_id = sample_id
  }

  output {
    File bam_error_rates = BamErrorRates.error_rates
  }
}

task BamErrorRates {
  input {
    File input_bam
    String sample_id
  }

  command <<<
    set -euo pipefail
    /usr/local/bin/bam_error_rates.py \
      -i ~{input_bam} \
      -o ~{sample_id}.bam.error_rates.tsv
  >>>

  output {
    File error_rates = "~{sample_id}.bam.error_rates.tsv"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/bam_error_rates"
    cpu: 1
    memory: "4G"
    disks: "local-disk 50 HDD"
  }
}