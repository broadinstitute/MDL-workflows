version 1.0


workflow bamToFastq {
    input {
        File inputBAM
        String sampleName
        Int maxRetries
     }

     call bamToFastq_task as b2fq {
         input:
           inputBAM=inputBAM,
           sampleName=sampleName,
           maxRetries=maxRetries
     }

     output {
         File fastq_file = b2fq.fastq_file
     }
}


task bamToFastq_task {
    input {
        File inputBAM
        String sampleName
        Int maxRetries
    }

    String output_fastq_filename = sampleName + ".fastq.gz"

    command <<<

        set -ex 
        set -o pipefail
        
        samtools fastq ~{inputBAM} | gzip -c > ~{output_fastq_filename}

        
    >>>

    output {
        File fastq_file = "~{output_fastq_filename}"
    }

    runtime {
        docker: "mgibio/samtools:1.16.1"
        cpu: 1
        memory: "4GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*2 + 5) + " SSD"
        preemptible: maxRetries
    }
}
