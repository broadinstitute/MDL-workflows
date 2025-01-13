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
        Int cpu = 8
        Int memoryGB = 32
        Int diskSizeGB
        Int preemptible_tries
    }

    String docker = "trinityctat/minimap2_lr:latest"

    String extra_arg = if allowSecondary then "" else "--secondary=no"
    String extra_arg2 = if keepUnmapped then "" else "--sam-hit-only"
    String extra_arg3 = if keepComments then "-y" else ""

    String extract_tags = if defined(tagsToExtract) && tagsToExtract != "" then "-T ~{tagsToExtract}" else ""

    command <<<
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

        fastq_name="temp.fastq"
        if [[ "~{inputExtension}" == "bam" ]] || [[ "$file_extension" == "ubam" ]]; then
            samtools fastq ~{extract_tags} ~{inputFile} > temp.fastq
        elif [[ "~{inputExtension}" == "fastq.gz" ]]; then
            mv ~{inputFile} temp.fastq.gz
            fastq_name="temp.fastq.gz"
        elif [[ "~{inputExtension}" == "fastq" ]]; then
            mv ~{inputFile} temp.fastq
        fi

        juncbed_arg=~{if defined(juncBED) then '"--junc-bed ${juncBED}"' else '""'}

        minimap2 ~{extra_arg2} ~{extra_arg3} -ax ${minimap2_preset} ~{customArguments} ${juncbed_arg} ~{extra_arg} -t ~{cpu} ~{referenceGenome} ${fastq_name} > temp.sam

        samtools sort -@ ~{cpu} temp.sam > ~{sampleName}.aligned.sorted.bam
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

workflow Minimap2_LR {
    meta {
        description: "Run Minimap2 from an unaligned BAM of PacBio long reads to generate an aligned sorted BAM and BAM index."
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
        Int cpu = 8
        Int memoryGB = 32
        Int ?diskSizeGB
        Int preemptible_tries = 3
    }


    String file_extension = basename(inputReads)
    if (sub(file_extension, "fastq$", "") != file_extension) {

        Int effective_diskSizeGB_fastq = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*4 + size(referenceGenome, "GB") + 20)])

        call Minimap2Task as minimap2_fastq {
            input:
                inputFile = inputReads,
                inputExtension = "fastq",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                tagsToExtract = tagsToExtract,
                keepComments = keepComments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                diskSizeGB = effective_diskSizeGB_fastq,
                cpu = cpu,
                memoryGB = memoryGB,
                preemptible_tries = preemptible_tries
        }
    }
    if (sub(file_extension, "fastq.gz$", "") != file_extension) {

        Int effective_diskSizeGB_fastqgz = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*20 + size(referenceGenome, "GB") + 20)])

        call Minimap2Task as minimap2_fastqgz {
            input:
                inputFile = inputReads,
                inputExtension = "fastq.gz",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                tagsToExtract = tagsToExtract,
                keepComments = keepComments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                diskSizeGB = effective_diskSizeGB_fastqgz,
                cpu = cpu,
                memoryGB = memoryGB,
                preemptible_tries = preemptible_tries
        }
    }
    if (sub(file_extension, "bam$", "") != file_extension) {
        
        Int effective_diskSizeGB_bam = select_first([diskSizeGB,  ceil(size(inputReads, "GB")*20 + size(referenceGenome, "GB") + 20)])
        
        call Minimap2Task as minimap2_ubam {
            input:
                inputFile = inputReads,
                inputExtension = "bam",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                tagsToExtract = tagsToExtract,
                keepComments = keepComments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                diskSizeGB = effective_diskSizeGB_bam,
                cpu = cpu,
                memoryGB = memoryGB,
                preemptible_tries = preemptible_tries
        }
    }

    output {
        File minimap2_bam = select_first([minimap2_fastq.minimap2_bam, minimap2_fastqgz.minimap2_bam, minimap2_ubam.minimap2_bam])
        File minimap2_bam_index = select_first([minimap2_fastq.minimap2_bam_index, minimap2_fastqgz.minimap2_bam_index, minimap2_ubam.minimap2_bam_index])
        File alignment_flagstat = select_first([minimap2_fastq.alignment_flagstat, minimap2_fastqgz.alignment_flagstat, minimap2_ubam.alignment_flagstat])
    }
}