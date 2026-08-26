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
        Int cpu = 30
        Int sortThreads = 2
        String sortMemory = "768M"
        Int memoryGB = 30
        Int? diskSizeGB
        Int preemptible_tries = 3
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/minimap2/minimap2:2.30-slim"

    # C3D only ships fixed vCPU tiers (4/8/16/30/60/90/180/360). Round the requested cpu up to
    # the smallest tier that covers it, and separately round memoryGB up to the smallest tier
    # whose highmem variant covers it; the larger of the two tiers is what actually gets
    # provisioned. Benchmarking (2026-08-25, see
    # minimap2_LR_n2d_vs_c3d_benchmark_20260825.md) showed c3d-highcpu-30 beats the old
    # n2d-highcpu-48 default at less than half the SPOT cost, so paying for a full tier instead
    # of chasing an exact custom match is the right tradeoff.
    Int cpu_tier = if cpu <= 4 then 4
        else if cpu <= 8 then 8
        else if cpu <= 16 then 16
        else if cpu <= 30 then 30
        else if cpu <= 60 then 60
        else if cpu <= 90 then 90
        else if cpu <= 180 then 180
        else 360
    Int mem_tier = if memoryGB <= 32 then 4
        else if memoryGB <= 64 then 8
        else if memoryGB <= 128 then 16
        else if memoryGB <= 240 then 30
        else if memoryGB <= 480 then 60
        else if memoryGB <= 720 then 90
        else if memoryGB <= 1440 then 180
        else 360
    Int effective_cpu = if cpu_tier >= mem_tier then cpu_tier else mem_tier
    # c3d-highcpu RAM per tier isn't a clean 2 GB/vCPU multiple at every size (e.g. 59 GB, not
    # 60, at the 30-vCPU tier), so the real values are listed explicitly.
    Int highcpu_ram = if effective_cpu == 4 then 8
        else if effective_cpu == 8 then 16
        else if effective_cpu == 16 then 32
        else if effective_cpu == 30 then 59
        else if effective_cpu == 60 then 118
        else if effective_cpu == 90 then 177
        else if effective_cpu == 180 then 354
        else 708
    String machine_type = if memoryGB <= highcpu_ram
        then "c3d-highcpu-~{effective_cpu}"
        else if memoryGB <= effective_cpu * 4
        then "c3d-standard-~{effective_cpu}"
        else "c3d-highmem-~{effective_cpu}"

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
        minimap2 ~{extra_arg2} ~{extra_arg3} -a ${preset_arg} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ~{effective_cpu} ~{referenceGenome} '~{sep="' '" inputFastqs}' \
            | samtools sort --no-PG --write-index -@ ~{sortThreads} -m ~{sortMemory} -O BAM \
                -o ~{sorted_bam_name}##idx##~{sorted_bam_index_name} -

        samtools flagstat ~{sorted_bam_name} > ~{sampleName}_alignment.flagstat.txt
    >>>

    output {
        File minimap2_bam = "~{sorted_bam_name}"
        File minimap2_bam_index = "~{sorted_bam_index_name}"
        File alignment_flagstat = "~{sampleName}_alignment.flagstat.txt"
    }

    # cpu/memory intentionally omitted: GCP Batch has a known bug where specifying both
    # predefinedMachineType and an explicit compute_resource (cpu_milli/memory_mib) can
    # spuriously reject an otherwise-valid combination ("machine_type ... cannot satisfy
    # compute_resource ..."), even when the values exactly match the machine type's real spec.
    # Let predefinedMachineType alone determine the shape.
    runtime {
        predefinedMachineType: "~{machine_type}"
        disks: "local-disk ~{effective_disk} SSD"
        docker: docker
        checkpointFiles: "monitoring.log"
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
        Int cpu = 30
        Int sortThreads = 2
        String sortMemory = "768M"
        Int memoryGB = 30
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
            sortThreads = sortThreads,
            sortMemory = sortMemory,
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
