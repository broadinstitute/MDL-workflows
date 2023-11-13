version 1.0

task Minimap2Task {
    input {
        File inputBAM
        File juncBED
        File referenceGenome
        String sampleName
        Boolean keepUnmapped = true
        Boolean allowSecondary = true
        Int cpu = 8
        Int memoryGB = 32
        Int diskSizeGB = 200
        String docker = "trinityctat/minimap2_lr:latest"
        Int preemptible_tries
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    String extra_arg = if allowSecondary then "" else "--secondary=no"
    String extra_arg2 = if keepUnmapped then "" else "--sam-hit-only"

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        minimap2 ~{extra_arg2} -ax splice:hq -uf --junc-bed ~{juncBED} ~{extra_arg} -t ~{cpu}  -G 1000 ~{referenceGenome} temp.fastq > temp.sam

        samtools sort -@ ~{cpu} temp.sam > ~{sampleName}.aligned.sorted.bam
        samtools index -@ ~{cpu} ~{sampleName}.aligned.sorted.bam
    >>>

    output {
        File minimap2_bam = "~{sampleName}.aligned.sorted.bam"
        File minimap2_bam_index = "~{sampleName}.aligned.sorted.bam.bai"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} SSD"
        docker: docker
        preemptible: preemptible_tries
    }
}

workflow Minimap2_LR {
    meta {
        description: "Run Minimap2 from an unaligned BAM of PacBio long reads to generate an aligned sorted BAM and BAM index."
    }

    input {
        File inputBAM
        File referenceGenome
        File juncBED
        String sampleName
        Boolean keepUnmapped = true
        Boolean allowSecondary = true
        Int preemptible_tries = 3
    }

    call Minimap2Task {
        input:
            inputBAM = inputBAM,
            juncBED = juncBED,
            referenceGenome = referenceGenome,
            sampleName = sampleName,
            keepUnmapped = keepUnmapped,
            allowSecondary = allowSecondary,
            preemptible_tries = preemptible_tries
    }

    output {
        File minimap2_bam = Minimap2Task.minimap2_bam
        File minimap2_bam_index = Minimap2Task.minimap2_bam_index
        File monitoringLog = Minimap2Task.monitoringLog
    }
}