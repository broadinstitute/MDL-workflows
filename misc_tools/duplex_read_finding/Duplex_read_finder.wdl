version 1.0

workflow Duplex_read_finder_wf {

    input {
        String sample_id
        File fastq
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/duplex-read-finding/duplex-read-finding:latest"
        Int preemptible = 0
    }
    
    call Duplex_read_finder_task as drf {

        input:
        sample_id=sample_id,
        fastq=fastq,
        preemptible=preemptible,
        docker=docker
    }

    output {
        File duplex_read_report = drf.duplex_read_report
    }
}



task Duplex_read_finder_task {

    input {
        String sample_id
        File fastq
        String docker
        Int preemptible = 0
    }

    String duplex_read_report_filename = "~{sample_id}.duplex_read_report.tsv.gz"
    Int disk_space_multiplier = 3
    Int disk_space = ceil(size(fastq, "GB")*disk_space_multiplier)


    command <<<

        set -ex
        set -o pipefail

        /usr/local/bin/find_duplex_reads.py --fastq ~{fastq} | gzip -c > ~{duplex_read_report_filename}

     >>>

     output {
         File duplex_read_report = "~{duplex_read_report_filename}"
     }


     runtime {
         docker: docker
         memory: "8GB"
         bootDiskSizeGb: 12
         disks: "local-disk ~{disk_space} HDD"
         cpu: 1
         preemptible: preemptible
     }
     
}

 
    
