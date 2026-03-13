version 1.0

import "minimap2_LR_fastq_list.wdl" as minimap2_input_list


workflow Minimap2_LR_fastq_list_from_set {
    meta {
        description: "Sample-set wrapper to run one Minimap2 alignment from multiple FASTQs."
    }

    input {
        Array[String] sample_names
        Array[Array[File]] fastq_pairs
        File referenceGenome
        File ?juncBED
        String readType
        String ?customArguments
        Boolean keepComments = true
        Boolean keepUnmapped = true
        Boolean allowSecondary = true
        Int cpu = 48
        Int memoryGB = 48
        Int preemptible_tries = 3
    }

    Array[Int] lane_indexes = range(length(sample_names))

    scatter (i in lane_indexes) {
        call minimap2_input_list.Minimap2MultiFastqTask as align {
            input:
                inputFastqs = fastq_pairs[i],
                referenceGenome = referenceGenome,
                juncBED = juncBED,
                sampleName = sample_names[i],
                readType = readType,
                customArguments = customArguments,
                keepComments = keepComments,
                keepUnmapped = keepUnmapped,
                allowSecondary = allowSecondary,
                cpu = cpu,
                memoryGB = memoryGB,
                preemptible_tries = preemptible_tries
        }
    }

    output {
        Array[File] minimap2_bams = align.minimap2_bam
        Array[File] minimap2_bam_indexes = align.minimap2_bam_index
        Array[File] alignment_flagstats = align.alignment_flagstat
    }
}
