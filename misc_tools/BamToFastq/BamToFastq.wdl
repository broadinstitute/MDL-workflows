version 1.0

workflow BamToFastq_wf {

    input {
        File input_bam
        String? tags_to_extract
        Int? exclude_flags
        Int compression_level = 6
        Int preemptible = 0

    }

    call BamToFastq_task {
        input:
          bam_file = input_bam,
          tags_to_extract = tags_to_extract,
          exclude_flags = exclude_flags,
          compression_level = compression_level,
          preemptible = preemptible
    }

    output {
        File fastq_gz = BamToFastq_task.fastq_gz
    }

}



task BamToFastq_task {
  input {
    File bam_file
    String? tags_to_extract
    # samtools -F: skip reads with any of these flag bits. Unset keeps the samtools default (0x900 = secondary + supplementary)
    Int? exclude_flags
    # samtools default is 1 (fast, ~50% larger output); 6 is near gzip -6 size at little time cost, 9 is very slow
    Int compression_level = 6
    Int preemptible = 0
  }

  String extract_tags = if defined(tags_to_extract) && tags_to_extract != "" then "-T ~{tags_to_extract}" else ""


  String exclude_arg = if defined(exclude_flags) then "-F ~{exclude_flags}" else ""

  Int disk_space_multiplier = 5
  Int disk_space = ceil(size(bam_file, "GB")*disk_space_multiplier)

      
  command {
    set -euo pipefail

    # Strip .bam extension and add .fastq.gz
    # -o/-0 to the same file captures both flagged-paired and unflagged reads;
    # samtools compresses the .gz output itself using the -@ thread pool (4 = vCPUs of c3d-highcpu-4)
    samtools fastq -@ 4 -c ~{compression_level} ~{extract_tags} ~{exclude_arg} -o ~{basename(bam_file, ".bam")}.fastq.gz -0 ~{basename(bam_file, ".bam")}.fastq.gz ~{bam_file}
    
  }

  output {
    File fastq_gz = "~{basename(bam_file, ".bam")}.fastq.gz"
  }

  runtime {
        predefinedMachineType: "c3d-highcpu-4"
        docker:"us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} SSD"
        preemptible: preemptible
    }
   
}