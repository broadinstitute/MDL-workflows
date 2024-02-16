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
    String juncbed_arg = if defined(juncBED) then "--juncBED ~{juncBED}" else ""


    command <<<
        bash ~{monitoringScript} > monitoring.log &

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
        elif [ "~{readType}" == "None" ]; then
            minimap2_preset=""
        else
            echo "Invalid readType: ~{readType}"
            exit 1
        fi

        fastq_name="temp.fastq"
        if [[ "~{inputExtension}" == "bam" ]] || [[ "$file_extension" == "ubam" ]]; then
            samtools fastq ~{inputFile} > temp.fastq
        elif [[ "~{inputExtension}" == "fastq.gz" ]]; then
            mv ~{inputFile} temp.fastq.gz
            fastq_name="temp.fastq.gz"
        elif [[ "~{inputExtension}" == "fastq" ]]; then
            mv ~{inputFile} temp.fastq
        fi

        minimap2 ~{extra_arg2} -ax ${minimap2_preset} ~{customArguments} ~{juncbed_arg} ~{extra_arg} -t ~{cpu} ~{referenceGenome} ${fastq_name} > temp.sam

        # minimap2 ~{extra_arg2} -ax splice:hq -uf --junc-bed ~{juncBED} ~{extra_arg} -t ~{cpu}  -G 1000 ~{referenceGenome} temp.fastq > temp.sam

        samtools sort -@ ~{cpu} temp.sam > ~{sampleName}.aligned.sorted.bam
        samtools index -@ ~{cpu} ~{sampleName}.aligned.sorted.bam

        samtools flagstat ~{sampleName}.aligned.sorted.bam > ~{sampleName}_alignment.flagstat.txt
    >>>

    output {
        File minimap2_bam = "~{sampleName}.aligned.sorted.bam"
        File minimap2_bam_index = "~{sampleName}.aligned.sorted.bam.bai"
        File alignment_flagstat = "~{sampleName}_alignment.flagstat.txt"
        File monitoringLog = "monitoring.log"
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
        Boolean keepUnmapped = true
        Boolean allowSecondary = false
        Int preemptible_tries = 3
    }


    String file_extension = basename(inputReads)
    if (sub(file_extension, "fastq$", "") != file_extension) {
        call Minimap2Task as minimap2_fastq {
            input:
                inputFile = inputReads,
                inputExtension = "fastq",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                preemptible_tries = preemptible_tries
        }
    }
    if (sub(file_extension, "fastq.gz$", "") != file_extension) {
        call Minimap2Task as minimap2_fastqgz {
            input:
                inputFile = inputReads,
                inputExtension = "fastq.gz",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                preemptible_tries = preemptible_tries
        }
    }
    if (sub(file_extension, "bam$", "") != file_extension) {
        call Minimap2Task as minimap2_ubam {
            input:
                inputFile = inputReads,
                inputExtension = "bam",
                juncBED = juncBED,
                referenceGenome = referenceGenome,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                preemptible_tries = preemptible_tries
        }
    }

    #File minimap2_bam = select_first(minimap2_fastq.minimap2_bam, minimap2_fastqgz.minimap2_bam, minimap2_ubam.minimap2_bam)
    #File minimap2_bam_index = select_first(minimap2_fastq.minimap2_bam_index, minimap2_fastqgz.minimap2_bam_index, minimap2_ubam.minimap2_bam_index)
    #File monitoringLog = select_first(minimap2_fastq.monitoringLog, minimap2_fastqgz.monitoringLog, minimap2_ubam.minimap2_bam_index)

    output {
        # File minimap2_bam = Minimap2Task.minimap2_bam
        File minimap2_bam = select_first([minimap2_fastq.minimap2_bam, minimap2_fastqgz.minimap2_bam, minimap2_ubam.minimap2_bam])
        # File minimap2_bam_index = Minimap2Task.minimap2_bam_index
        File minimap2_bam_index = select_first([minimap2_fastq.minimap2_bam_index, minimap2_fastqgz.minimap2_bam_index, minimap2_ubam.minimap2_bam_index])
        File alignment_flagstat = select_first([minimap2_fastq.alignment_flagstat, minimap2_fastqgz.alignment_flagstat, minimap2_ubam.alignment_flagstat])
        # File monitoringLog = Minimap2Task.monitoringLog
        File monitoringLog = select_first([minimap2_fastq.monitoringLog, minimap2_fastqgz.monitoringLog, minimap2_ubam.minimap2_bam_index])
    }
}