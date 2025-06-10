version 1.0

import "minimap2_wrapper.wdl" as minimap2_wrapper


task splitReadsTask {
    input {
        File inputFile
        String sampleName  # sampleName
        Int reads_per_split
        Int preemptible_tries = 3
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/split-reads-in-chunks/split-reads-in-chunks:latest"
    Int diskSizeGB = ceil(size(inputFile, "GB")*2.2) + 20

    command <<<

        # adding _splitReads to the output prefix so we can make the output glob more precise and avoid globbing the input back
        python /scripts/split_reads_in_chunks.py \
            --input_file ~{inputFile} \
            --output_prefix ~{sampleName}_splitReads \
            --chunk_size ~{reads_per_split}

    >>>

    output {
        Array[File] split_inputs = glob("~{sampleName}_splitReads_*")
    }

    runtime {
        cpu: 4
        memory: "4 GB"
        disks: "local-disk ~{diskSizeGB} SSD"
        docker: docker
        preemptible: preemptible_tries
    }
}


task mergeBAMs {
    input {
        Array[File] bams_to_merge
        String sampleName
        Int preemptible_tries
    }

    Int diskSizeGB = ceil(size(bams_to_merge, "GB") * 3) + 20

    command <<<

        samtools merge --threads 4 -o ~{sampleName}.bam '~{sep="' '" bams_to_merge}'
        samtools index  ~{sampleName}.bam
        samtools flagstats  ~{sampleName}.bam > ~{sampleName}.flagstat.txt

    >>>

    output {
        File merged_bam = "~{sampleName}.bam"
        File merged_bam_index = "~{sampleName}.bam.bai"
        File alignment_flagstat = "~{sampleName}.flagstat.txt"
    }

    runtime {
        cpu: 4
        memory: "4 GB"
        disks: "local-disk ~{diskSizeGB} SSD"
        docker: "mgibio/samtools:v1.21-noble"
        preemptible: preemptible_tries
    }
}


workflow Minimap2_LR {
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
        Int? reads_per_shard
        Int cpu = 8
        Int memoryGB = 32
        Int ?diskSizeGB
        Int preemptible_tries = 3
    }


    if (defined(reads_per_shard)) {
        call splitReadsTask {
            input:
                inputFile = inputReads,
                sampleName = sampleName,
                reads_per_split = select_first([reads_per_shard]),
                preemptible_tries = preemptible_tries
        }


        scatter(split_input in splitReadsTask.split_inputs) {
                call minimap2_wrapper.Minimap2_wrapper as minimap2_shard {
                    input:
                        inputReads = split_input,
                        referenceGenome = referenceGenome,
                        juncBED = juncBED,
                        sampleName = basename(split_input),
                        readType = readType,
                        customArguments = customArguments,
                        tagsToExtract = tagsToExtract,
                        keepComments = keepComments,
                        keepUnmapped = keepUnmapped,
                        allowSecondary = allowSecondary,
                        cpu = cpu,
                        memoryGB = memoryGB,
                        preemptible_tries = preemptible_tries
                }
        }

        call mergeBAMs {
            input:
                bams_to_merge = minimap2_shard.minimap2_bam,
                sampleName = sampleName,
                preemptible_tries = preemptible_tries
        }
    }


    if (!defined(reads_per_shard)) {
        call minimap2_wrapper.Minimap2_wrapper as minimap2_whole {
            input:
                inputReads = inputReads,
                referenceGenome = referenceGenome,
                juncBED = juncBED,
                sampleName = sampleName,
                readType = readType,
                customArguments = customArguments,
                tagsToExtract = tagsToExtract,
                keepComments = keepComments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                cpu = cpu,
                memoryGB = memoryGB,
                preemptible_tries = preemptible_tries
        }
    }

 

    output {
        File minimap2_bam = select_first([mergeBAMs.merged_bam, minimap2_whole.minimap2_bam])
        File minimap2_bam_index = select_first([mergeBAMs.merged_bam_index, minimap2_whole.minimap2_bam_index])
        File alignment_flagstat = select_first([mergeBAMs.alignment_flagstat, minimap2_whole.alignment_flagstat])
    }
}