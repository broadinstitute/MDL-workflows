version 1.0

workflow CUDLL_scattered {
    input {
        Array[File] input_bams
        Array[File]? input_bais
        File? reference_fasta
        File? reference_fai

        String sample_name
        String barcode_tag = "CB"
        String umi_tag = "UB"
        Float identity = 0.95
        String? tags
        Boolean no_consensus = false
        Boolean emit_supplementary_alignments = true
        Boolean emit_consensus_sorted = false

        Int? cpu
        Int? memory_gb

        String docker_image_cudll
    }

    Boolean need_index = !defined(input_bais) || length(select_first([input_bais])) != length(input_bams)

    scatter (bam_idx in range(length(input_bams))) {
        File input_bam = input_bams[bam_idx]
        String shard_prefix = basename(input_bam, ".bam")

        if (need_index) {
            call CreateIndex { input: input_bam = input_bam, docker_image = docker_image_cudll }
        }

        # Use provided index if available and length matches, otherwise use generated index
        File bam_index = select_first([
            if defined(input_bais) && length(select_first([input_bais])) == length(input_bams)
                then select_first([input_bais])[bam_idx]
                else CreateIndex.bam_index
        ])

        call LocalOverlap {
            input:
                input_bam = input_bam,
                input_bai = bam_index,
                reference_fasta = reference_fasta,
                reference_fai = reference_fai,
                output_prefix = shard_prefix,
                barcode_tag = barcode_tag,
                umi_tag = umi_tag,
                tags = tags,
                no_consensus = no_consensus,
                emit_supplementary_alignments = emit_supplementary_alignments,
                emit_consensus_sorted = emit_consensus_sorted,
                cpu = cpu,
                memory_gb = memory_gb,
                docker_image = docker_image_cudll
        }

        call CrossLocus {
            input:
                consensus_bam = LocalOverlap.consensus_bam,
                output_prefix = shard_prefix,
                barcode_tag = barcode_tag,
                umi_tag = umi_tag,
                identity = identity,
                cpu = cpu,
                memory_gb = memory_gb,
                docker_image = docker_image_cudll
        }
    }

    call MergeFinalBams {
        input:
            bams = CrossLocus.final_bam,
            output_name = sample_name + ".merged.bam"
    }

    # Merge supplementary alignment BAMs if they exist (they are not guaranteed to be sorted)
    if (emit_supplementary_alignments) {
        Array[File] supplementary_bams_filtered = select_all(LocalOverlap.supplementary_alignments_bam)
        if (length(supplementary_bams_filtered) > 0) {
            call MergeUnsortedBams as MergeSupplementaryBams {
                input:
                    bams = supplementary_bams_filtered,
                    output_name = sample_name + ".supplementary_alignments.merged.bam"
            }
        }
    }

    output {
        File   merged_bam                       = MergeFinalBams.merged_bam
        File   merged_bai                       = MergeFinalBams.merged_bai
        File?  merged_supplementary_bam         = MergeSupplementaryBams.merged_bam
        File?  merged_supplementary_bai         = MergeSupplementaryBams.merged_bai
        Array[File]  shard_final_bams           = CrossLocus.final_bam
        Array[File]  shard_final_bais           = CrossLocus.final_bai
        Array[File?] shard_supplementary_bams   = LocalOverlap.supplementary_alignments_bam
        Array[File?] shard_consensus_sorted_bam = LocalOverlap.consensus_sorted_bam
        Array[File?] shard_consensus_sorted_bai = LocalOverlap.consensus_sorted_bai
    }
}

task CreateIndex {
    input {
        File input_bam
        String docker_image
    }

    Int disk_gb = ceil(size(input_bam, "GB")) + 10

    command <<<
        set -euo pipefail
        samtools index -@ 2 "~{input_bam}"
    >>>

    output {
        File bam_index = "~{input_bam}.bai"
    }

    runtime {
        cpu: 2
        memory: "2 GB"
        docker: docker_image
        disks: "local-disk ~{disk_gb} SSD"
        predefinedMachineType: "n1-highcpu-2"
        preemptible: 3
    }
}

task LocalOverlap {
    input {
        File input_bam
        File? input_bai
        File? reference_fasta
        File? reference_fai

        String output_prefix
        String barcode_tag
        String umi_tag
        String? tags
        Boolean no_consensus
        Boolean emit_supplementary_alignments
        Boolean emit_consensus_sorted

        Int? cpu
        Int? memory_gb

        String docker_image
    }

    Int task_cpu = select_first([cpu, 32])
    Int task_memory_gb = select_first([memory_gb, 32])
    String? machine_type = if defined(cpu) || defined(memory_gb) then "n2d-custom-~{task_cpu}-~{task_memory_gb * 1024}" else "n2d-highcpu-32"

    String tags_arg = if defined(tags) then "--tags " + tags else ""
    String supplementary_alignments_bam_path = output_prefix + ".supplementary_alignments.bam"
    String supplementary_alignments_arg = if emit_supplementary_alignments then "--sa-read-bam " + supplementary_alignments_bam_path else ""
    Int disk_gb = ceil(size(input_bam, "GB") * (if emit_consensus_sorted then 3 else 2)) + 20

    command <<<
        set -euo pipefail

        cudll_local_overlap \
            -i "~{input_bam}" \
            -o "~{output_prefix}.consensus.bam" \
            ~{if defined(reference_fasta) then "-r \"" + select_first([reference_fasta]) + "\"" else ""} \
            ~{tags_arg} \
            ~{supplementary_alignments_arg} \
            ~{if no_consensus then "--no-consensus" else ""} \
            -t ~{task_cpu} \
            --barcode ~{barcode_tag} \
            --umi ~{umi_tag} \
            --umi-hamming-only

        if [ "~{emit_consensus_sorted}" = "true" ]; then
            samtools sort -@ ~{task_cpu} \
            -o "~{output_prefix}.consensus.sorted.bam" \
            "~{output_prefix}.consensus.bam"
            samtools index -@ ~{task_cpu} "~{output_prefix}.consensus.sorted.bam"
        fi
    >>>

    output {
        File consensus_bam = "~{output_prefix}.consensus.bam"
        File? supplementary_alignments_bam = supplementary_alignments_bam_path
        File? consensus_sorted_bam = "~{output_prefix}.consensus.sorted.bam"
        File? consensus_sorted_bai = "~{output_prefix}.consensus.sorted.bam.bai"
    }

    runtime {
        cpu: task_cpu
        memory: "~{task_memory_gb} GB"
        docker: docker_image
        disks: "local-disk ~{disk_gb} SSD"
        predefinedMachineType: machine_type
        preemptible: 3
    }
}

task CrossLocus {
    input {
        File consensus_bam
        String output_prefix
        String barcode_tag
        String umi_tag
        Float identity

        Int? cpu
        Int? memory_gb

        String docker_image
    }

    Int task_cpu = select_first([cpu, 32])
    Int task_memory_gb = select_first([memory_gb, 32])
    String? machine_type = if defined(cpu) || defined(memory_gb) then "n2d-custom-~{task_cpu}-~{task_memory_gb * 1024}" else "n2d-highcpu-32"

    Int disk_gb = ceil(size(consensus_bam, "GB") * 7) + 20

    command <<<
        set -euo pipefail

        samtools sort -@ ~{task_cpu} -t ~{barcode_tag} \
            -o "~{output_prefix}.consensus.~{barcode_tag}_sorted.bam" \
            "~{consensus_bam}"

        cudll_cross_locus \
            -i "~{output_prefix}.consensus.~{barcode_tag}_sorted.bam" \
            -o "~{output_prefix}.consensus.homology_dedup.~{barcode_tag}_sorted.bam" \
            -t ~{task_cpu} \
            --barcode ~{barcode_tag} \
            --umi ~{umi_tag} \
            --identity ~{identity} \
            --umi-hamming-only \
            --rank-by-aligned-bases

        samtools sort -@ ~{task_cpu} \
            -o "~{output_prefix}.consensus.homology_dedup.sorted.bam" \
            "~{output_prefix}.consensus.homology_dedup.~{barcode_tag}_sorted.bam"

        samtools index -@ ~{task_cpu} "~{output_prefix}.consensus.homology_dedup.sorted.bam"
    >>>

    output {
        File final_bam = "~{output_prefix}.consensus.homology_dedup.sorted.bam"
        File final_bai = "~{output_prefix}.consensus.homology_dedup.sorted.bam.bai"
    }

    runtime {
        cpu: task_cpu
        memory: "~{task_memory_gb} GB"
        docker: docker_image
        disks: "local-disk ~{disk_gb} SSD"
        predefinedMachineType: machine_type
        preemptible: 3
    }
}

task MergeFinalBams {
    input {
        Array[File] bams
        String output_name
    }

    Int diskGB = ceil(size(bams, "GB") * 2.5 + 20)

    command <<<
        set -euo pipefail

        samtools merge -@ 2 -o ~{output_name} ~{sep=' ' bams}
        samtools index -@ 2 ~{output_name}
    >>>

    output {
        File merged_bam = "~{output_name}"
        File merged_bai = "~{output_name}.bai"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        cpu: 2
        memory: "2 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 2
        predefinedMachineType: "n2d-highcpu-2"
    }
}

task MergeUnsortedBams {
    input {
        Array[File] bams
        String output_name
    }

    Int diskGB = ceil(size(bams, "GB") * 3 + 20)

    command <<<
        set -euo pipefail

        # Sort each BAM by coordinate before merging
        for bam in ~{sep=' ' bams}; do
            samtools sort -@ 2 -o "sorted_$(basename "$bam")" "$bam"
        done

        # Merge sorted BAMs and index
        samtools merge -@ 2 -o ~{output_name} sorted_*.bam
        samtools index -@ 2 ~{output_name}
    >>>

    output {
        File merged_bam = "~{output_name}"
        File merged_bai = "~{output_name}.bai"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        cpu: 2
        memory: "2 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 2
        predefinedMachineType: "n2d-highcpu-2"
    }
}
