version 1.0

workflow BamToFastq_wf {

    input {
        File input_bam
        Int preemptible = 0

    }

    call BamToFastq_task {
        input:
          bam_file = input_bam,
          preemptible = preemptible
    }

    output {
        File fastq_gz = BamToFastq_task.fastq_gz
    }

}



task BamToFastq_task {
  input {
    File bam_file
    Int preemptible = 0
  }


  Int disk_space_multiplier = 5
  Int disk_space = ceil(size(bam_file, "GB")*disk_space_multiplier)

      
  command {
    set -euo pipefail

    # Strip .bam extension and add .fastq.gz
    samtools fastq ~{bam_file} | gzip > ~{basename(bam_file, ".bam")}.fastq.gz
    
  }

  output {
    File fastq_gz = "~{basename(bam_file, ".bam")}.fastq.gz"
  }

  runtime {
        docker:"us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        memory: "8GB"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
   
}