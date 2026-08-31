version 1.0

workflow MergeBams {
    input {
        Array[File]+ input_bams
        String output_name
    }

    call Merge {
        input:
            bams = input_bams,
            output_name = output_name
    }

    output {
        File merged_bam = Merge.merged_bam
        File merged_bai = Merge.merged_bai
    }
}

task Merge {
    input {
        Array[File]+ bams
        String output_name
    }

    Int disk_gb = ceil(size(bams, "GB") * 2.5) + 20

    command <<<
        set -euo pipefail
        samtools merge --no-PG --write-index -p -@ 4 \
            -o ~{output_name}##idx##~{output_name}.bai \
            ~{sep=' ' bams}
    >>>

    output {
        File merged_bam = "~{output_name}"
        File merged_bai = "~{output_name}.bai"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: 2
        predefinedMachineType: "c3d-highcpu-4"
    }
}
