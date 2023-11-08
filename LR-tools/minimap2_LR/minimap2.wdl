version 1.0

task Minimap2Task {
    input {
        File inputBAM
        File juncBED
        File referenceGenome
        String sampleName
        Int cpu = 8
        Int memoryGB = 32
        Int diskSizeGB = 200
        String docker = "trinityctat/minimap2_lr:latest"
        Int preemptible_tries
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools fastq ~{inputBAM} > temp.fastq

        minimap2 --sam-hit-only -ax splice:hq -uf --junc-bed ~{juncBED} --secondary=no -t ~{cpu}  -G 1000 ~{referenceGenome} temp.fastq > temp.sam

        samtools sort temp.sam > ~{sampleName}.aligned.sorted.bam
        samtools index ~{sampleName}.aligned.sorted.bam
    >>>

    output {
        File minimap2_bam = "~{sampleName}.aligned.sorted.bam"
        File minimap2_bam_index = "~{sampleName}.aligned.sorted.bam.bai"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        preemptible: preemptible_tries
    }
}

workflow Minimap2_LR {
    input {
        File inputBAM
        File referenceGenome
        File referenceGenomeIndex
        File juncBED
        String sampleName
        Int preemptible_tries = 3
    }

    call Minimap2Task {
        input:
            inputBAM = inputBAM,
            juncBED = juncBED,
            referenceGenome = referenceGenome,
            sampleName = sampleName,
            preemptible_tries = preemptible_tries
    }

    output {
        File minimap2_bam = Minimap2Task.minimap2_bam
        File minimap2_bam_index = Minimap2Task.minimap2_bam_index
        File monitoringLog = Minimap2Task.monitoringLog
    }
}