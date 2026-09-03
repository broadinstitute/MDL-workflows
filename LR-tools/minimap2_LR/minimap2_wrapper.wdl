version 1.0


task Minimap2Task {
    input {
        File inputFile
        String inputExtension
        File ?juncBED
        File referenceGenome
        String sampleName
        String readType
        String ?customArguments
        String ?tagsToExtract
        Boolean keepComments = true
        Boolean keepUnmapped = true
        Boolean allowSecondary = true
        Boolean writeIndex = true
        Int cpu = 30
        Int sortThreads = 2
        String sortMemory = "768M"
        Int memoryGB = 59
        Int diskSizeGB
        Int preemptible_tries
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/minimap2/minimap2:latest"

    # C3D only ships fixed vCPU tiers (4/8/16/30/60/90/180/360). Round the requested cpu up to
    # the smallest tier that covers it, and separately round memoryGB up to the smallest tier
    # whose highmem variant covers it; the larger of the two tiers is what actually gets
    # provisioned. Benchmarking (2026-08-25, see minimap2_LR_n2d_vs_c3d.md) showed
    # c3d-highcpu-30 beats the old n2d-highcpu-48 default at less than half the SPOT cost, so
    # paying for a full tier instead of chasing an exact custom match is the right tradeoff.
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

    String extra_arg = if allowSecondary then "" else "--secondary=no"
    String extra_arg2 = if keepUnmapped then "" else "--sam-hit-only"
    String extra_arg3 = if keepComments then "-y" else ""

    String extract_tags = if defined(tagsToExtract) && tagsToExtract != "" then "-T ~{tagsToExtract}" else ""
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

        sort_args=(--no-PG -@ ~{sortThreads} -m ~{sortMemory} -O BAM)
        if [ "~{writeIndex}" == "true" ]; then
            sort_args+=(--write-index)
            sort_out="~{sorted_bam_name}##idx##~{sorted_bam_index_name}"
        else
            sort_out="~{sorted_bam_name}"
        fi

        if [[ "~{inputExtension}" == "bam" ]]; then
            samtools fastq ~{extract_tags} ~{inputFile} \
                | minimap2 ~{extra_arg2} ~{extra_arg3} -ax ${minimap2_preset} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ~{effective_cpu} ~{referenceGenome} - \
                | samtools sort "${sort_args[@]}" -o "$sort_out" -
        elif [[ "~{inputExtension}" == "fastq.zst" ]]; then
            zstd -d -c ~{inputFile} \
                | minimap2 ~{extra_arg2} ~{extra_arg3} -ax ${minimap2_preset} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ~{effective_cpu} ~{referenceGenome} - \
                | samtools sort "${sort_args[@]}" -o "$sort_out" -
        elif [[ "~{inputExtension}" == "fastq.gz" ]] || [[ "~{inputExtension}" == "fastq" ]]; then
            minimap2 ~{extra_arg2} ~{extra_arg3} -ax ${minimap2_preset} ~{custom_args} ~{if defined(juncBED) then "--junc-bed " + juncBED else ""} ~{extra_arg} -t ~{effective_cpu} ~{referenceGenome} ~{inputFile} \
                | samtools sort "${sort_args[@]}" -o "$sort_out" -
        else
            echo "Unsupported inputExtension: ~{inputExtension}"
            exit 1
        fi

        samtools flagstat ~{sorted_bam_name} > ~{sampleName}_alignment.flagstat.txt
    >>>

    output {
        File minimap2_bam = "~{sorted_bam_name}"
        File? minimap2_bam_index = "~{sorted_bam_index_name}"
        File alignment_flagstat = "~{sampleName}_alignment.flagstat.txt"
    }

    runtime {
        predefinedMachineType: "~{machine_type}"
        disks: "local-disk ~{diskSizeGB} SSD"
        docker: docker
        preemptible: preemptible_tries
    }
}



workflow Minimap2_wrapper {
    meta {
        description: "Run Minimap2 from an (unaligned) BAM or FASTQ of single end long reads to generate an aligned sorted BAM and BAM index."
    }

    input {
        File inputReads
        File referenceGenome
        File ?juncBED
        String sampleName
        String readType
        String ?customArguments
        String ?tagsToExtract
        Boolean keepComments = true
        Boolean keepUnmapped = true
        Boolean allowSecondary = false
        Boolean writeIndex = true
        Int cpu = 30
        Int sortThreads = 2
        String sortMemory = "768M"
        Int memoryGB = 59
        Int ?diskSizeGB
        Int preemptible_tries = 3
    }

    String file_name = basename(inputReads)

    if (sub(file_name, "fastq$", "") != file_name) {

        Int effective_diskSizeGB_fastq = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*4 + size(referenceGenome, "GB") + 20)])
        String inputExtension_fastq = "fastq"

    }
    if (sub(file_name, "fastq.gz$", "") != file_name) {

        Int effective_diskSizeGB_fastqgz = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*20 + size(referenceGenome, "GB") + 20)])
        String inputExtension_fastqgz = "fastq.gz"

    }
    if (sub(file_name, "fastq.zst$", "") != file_name) {

        Int effective_diskSizeGB_fastqzstd = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*25 + size(referenceGenome, "GB") + 20)])
        String inputExtension_fastqzstd = "fastq.zst"

    }
    if (sub(file_name, "bam$", "") != file_name) {
        
        Int effective_diskSizeGB_bam = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*20 + size(referenceGenome, "GB") + 20)])
        String inputExtension_bam = "bam"

    }

    Int effective_diskSizeGB = select_first([effective_diskSizeGB_fastq, effective_diskSizeGB_fastqgz, effective_diskSizeGB_fastqzstd, effective_diskSizeGB_bam])
    String inputExtension =  select_first([inputExtension_fastq, inputExtension_fastqgz, inputExtension_fastqzstd, inputExtension_bam])

    call Minimap2Task as minimap2_run {
        input:
            inputFile = inputReads,
            inputExtension = inputExtension,
            juncBED = juncBED,
            referenceGenome = referenceGenome,
            sampleName = sampleName,
            readType = readType,
            customArguments = customArguments,
            tagsToExtract = tagsToExtract,
            keepComments = keepComments,
            keepUnmapped = keepUnmapped,
            allowSecondary = allowSecondary,
            writeIndex = writeIndex,
            diskSizeGB = effective_diskSizeGB,
            cpu = cpu,
            sortThreads = sortThreads,
            sortMemory = sortMemory,
            memoryGB = memoryGB,
            preemptible_tries = preemptible_tries
    }

    output {
        File minimap2_bam = minimap2_run.minimap2_bam
        File? minimap2_bam_index = minimap2_run.minimap2_bam_index
        File alignment_flagstat = minimap2_run.alignment_flagstat
    }
}
