version 1.0

workflow Flagstat_wf {

    input {
        String sample_id
        File input_bam
        Int preemptible = 0

    }

    call Flagstat_task {
        input:
          sample_id = sample_id,
          input_bam = input_bam,
          preemptible = preemptible
    }

    output {
        File flagstat_file = Flagstat_task.flagstat_file
        String total_reads = Flagstat_task.total_reads
    }

}

task Flagstat_task {

    input {
        String sample_id
        File input_bam
        Int preemptible = 0
    }

    String stats_filename = "~{sample_id}.flagstat.txt"
    Int disk_space_multiplier = 2
    Int disk_space = ceil(size(input_bam, "GB")*disk_space_multiplier)
    
    command <<<
        set -ex

        samtools flagstat ~{input_bam} > ~{stats_filename}

        # Extract total reads count from primary alignments line (QC-passed + QC-failed)
        sed -n '2p' ~{stats_filename} | awk '{print $1 + $3}' > total_reads.txt

    >>>

    output {
        File flagstat_file = "~{stats_filename}"
        String total_reads = read_string("total_reads.txt")
    }

   runtime {
        docker:"quay.io/biocontainers/samtools:1.21--h96c455f_1"
        memory: "8GB"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
    
    
}
