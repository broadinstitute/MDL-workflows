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
        Int cpu = 8
        Int memoryGB = 32
        Int diskSizeGB
        Int preemptible_tries = 3
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/minimap2/minimap2:2.30-slim"

    String extra_arg = if allowSecondary then "" else "--secondary=no"
    String extra_arg2 = if keepUnmapped then "" else "--sam-hit-only"
    String extra_arg3 = if keepComments then "-y" else ""
    String custom_args = select_first([customArguments, ""])

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

        sort_threads=$(( ~{cpu} / 4 ))
        if [ "${sort_threads}" -lt 1 ]; then
            sort_threads=1
        fi
        if [ "${sort_threads}" -ge ~{cpu} ]; then
            sort_threads=$(( ~{cpu} - 1 ))
        fi
        if [ "${sort_threads}" -lt 0 ]; then
            sort_threads=0
        fi
        minimap2_threads=$(( ~{cpu} - sort_threads ))
        if [ "${minimap2_threads}" -lt 1 ]; then
            minimap2_threads=1
        fi

        minimap2 ~{extra_arg2} ~{extra_arg3} -a ${preset_arg} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ${minimap2_threads} ~{referenceGenome} '~{sep="' '" inputFastqs}' \
            | samtools sort -@ ${sort_threads} -O BAM -o ~{sampleName}.aligned.sorted.bam -

        samtools index -@ ~{cpu} ~{sampleName}.aligned.sorted.bam
        samtools flagstat ~{sampleName}.aligned.sorted.bam > ~{sampleName}_alignment.flagstat.txt
    >>>

    output {
        File minimap2_bam = "~{sampleName}.aligned.sorted.bam"
        File minimap2_bam_index = "~{sampleName}.aligned.sorted.bam.bai"
        File alignment_flagstat = "~{sampleName}_alignment.flagstat.txt"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GB"
        disks: "local-disk ~{diskSizeGB} SSD"
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
        Int cpu = 8
        Int memoryGB = 32
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
