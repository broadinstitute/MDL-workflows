version 1.0


task Minimap2MultiFastqTask {
    input {
        Array[File] inputFastqs
        File referenceGenome
        File ?juncBED
        String sampleName
        String readType
        String ?customArguments
        Boolean keepComments = true
        Boolean keepUnmapped = true
        Boolean allowSecondary = true
        Int cpu = 48
        Int memoryGB = 48
        Int? diskSizeGB
        Int preemptible_tries = 3
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/minimap2/minimap2:2.30-slim"
    Boolean use_predefined_machine_type = cpu == 2 || cpu == 4 || cpu == 8 || cpu == 16 || cpu == 32 || cpu == 48 || cpu == 64 || cpu == 80 || cpu == 96
    String machine_type = if (use_predefined_machine_type && memoryGB == cpu * 4)
        then "n2d-standard-~{cpu}"
        else if (use_predefined_machine_type && memoryGB == cpu * 8)
        then "n2d-highmem-~{cpu}"
        else if (use_predefined_machine_type && memoryGB == cpu)
        then "n2d-highcpu-~{cpu}"
        else "n2d-custom-~{cpu}-~{memoryGB * 1024}"

    Int effective_disk = select_first([diskSizeGB, ceil(size(inputFastqs, "GB") * 4 + size(referenceGenome, "GB") + 20)])
    String extra_arg = if allowSecondary then "" else "--secondary=no"
    String extra_arg2 = if keepUnmapped then "" else "--sam-hit-only"
    String extra_arg3 = if keepComments then "-y" else ""
    String custom_args = select_first([customArguments, ""])
    String sorted_bam_name = "~{sampleName}.aligned.sorted.bam"
    String sorted_bam_index_name = "~{sampleName}.aligned.sorted.bam.bai"

    command <<<
        set -euo pipefail

        minimap2_preset=""

        if [ "~{readType}" == "PacBioCLR" ]; then
            minimap2_preset="map-pb"
        elif [ "~{readType}" == "ONTGenomic" ]; then
            minimap2_preset="map-ont"
        elif [ "~{readType}" == "PacBioHiFi" ]; then
            minimap2_preset="map-hifi"
        elif [ "~{readType}" == "SplicedLongReads" ]; then
            minimap2_preset="splice"
        elif [ "~{readType}" == "ONTDirectRNA" ]; then
            minimap2_preset="splice -uf -k14"
        elif [ "~{readType}" == "PacBioIsoSeq" ]; then
            minimap2_preset="splice:hq -uf"
        elif [ "~{readType}" == "ONTGenomicQ20" ]; then
            minimap2_preset="lr:hq"
        elif [ "~{readType}" == "None" ]; then
            minimap2_preset=""
        else
            echo "Invalid readType: ~{readType}"
            exit 1
        fi

        if [ -n "${minimap2_preset}" ]; then
            preset_arg="-x ${minimap2_preset}"
        else
            preset_arg=""
        fi
        minimap2 ~{extra_arg2} ~{extra_arg3} -a ${preset_arg} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ~{cpu} ~{referenceGenome} '~{sep="' '" inputFastqs}' \
            | samtools sort --no-PG --write-index -@ ~{cpu} -O BAM \
                -o ~{sorted_bam_name}##idx##~{sorted_bam_index_name} -

        samtools flagstat ~{sorted_bam_name} > ~{sampleName}_alignment.flagstat.txt
    >>>

    output {
        File minimap2_bam = "~{sorted_bam_name}"
        File minimap2_bam_index = "~{sorted_bam_index_name}"
        File alignment_flagstat = "~{sampleName}_alignment.flagstat.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GB"
        predefinedMachineType: "~{machine_type}"
        disks: "local-disk ~{effective_disk} SSD"
        docker: docker
        preemptible: preemptible_tries
    }
}


workflow Minimap2_LR_fastq_list {
    meta {
        description: "Run one Minimap2 call against an ordered list of FASTQ files to produce a single aligned sorted BAM and BAM index."
    }

    input {
        Array[File] inputFastqs
        File referenceGenome
        File ?juncBED
        String sampleName
        String readType
        String ?customArguments
        Boolean keepComments = true
        Boolean keepUnmapped = true
        Boolean allowSecondary = false
        Int cpu = 48
        Int memoryGB = 48
        Int ?diskSizeGB
        Int preemptible_tries = 3
    }

    Int effective_diskSizeGB = select_first([diskSizeGB, ceil(size(inputFastqs, "GB") * 4 + size(referenceGenome, "GB") + 20)])

    call Minimap2MultiFastqTask as minimap2_run {
        input:
            inputFastqs = inputFastqs,
            referenceGenome = referenceGenome,
            juncBED = juncBED,
            sampleName = sampleName,
            readType = readType,
            customArguments = customArguments,
            keepComments = keepComments,
            keepUnmapped = keepUnmapped,
            allowSecondary = allowSecondary,
            cpu = cpu,
            memoryGB = memoryGB,
            diskSizeGB = effective_diskSizeGB,
            preemptible_tries = preemptible_tries
    }

    output {
        File minimap2_bam = minimap2_run.minimap2_bam
        File minimap2_bam_index = minimap2_run.minimap2_bam_index
        File alignment_flagstat = minimap2_run.alignment_flagstat
    }
}
