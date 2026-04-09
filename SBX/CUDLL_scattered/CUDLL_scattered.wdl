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
        Boolean prune_pg_header_merge_final_bams = false
        String mitochondrial_contig_name = "chrM"

        Int? cpu
        Int? memory_gb
        Int? local_overlap_mito_cpu
        Int? local_overlap_mito_memory_gb

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

        call LocalOverlap as LocalOverlapNonMito {
            input:
                input_bam = input_bam,
                input_bai = bam_index,
                reference_fasta = reference_fasta,
                reference_fai = reference_fai,
                output_prefix = shard_prefix + ".non_mito",
                barcode_tag = barcode_tag,
                umi_tag = umi_tag,
                tags = tags,
                no_consensus = no_consensus,
                emit_supplementary_alignments = emit_supplementary_alignments,
                emit_consensus_sorted = emit_consensus_sorted,
                mitochondrial_only = false,
                mitochondrial_contig_name = mitochondrial_contig_name,
                cpu = cpu,
                memory_gb = memory_gb,
                docker_image = docker_image_cudll
        }

        call LocalOverlap as LocalOverlapMito {
            input:
                input_bam = input_bam,
                input_bai = bam_index,
                reference_fasta = reference_fasta,
                reference_fai = reference_fai,
                output_prefix = shard_prefix + ".mito",
                barcode_tag = barcode_tag,
                umi_tag = umi_tag,
                tags = tags,
                no_consensus = no_consensus,
                emit_supplementary_alignments = emit_supplementary_alignments,
                emit_consensus_sorted = emit_consensus_sorted,
                mitochondrial_only = true,
                mitochondrial_contig_name = mitochondrial_contig_name,
                cpu = select_first([local_overlap_mito_cpu, 32]),
                memory_gb = if defined(local_overlap_mito_memory_gb) then local_overlap_mito_memory_gb else if defined(local_overlap_mito_cpu) then select_first([local_overlap_mito_cpu]) * 4 else if defined(memory_gb) then memory_gb * 4 else 128,
                docker_image = docker_image_cudll
        }

        call MergeTagSortedBams as MergeShardConsensusBams {
            input:
                bams = [LocalOverlapNonMito.consensus_bam, LocalOverlapMito.consensus_bam],
                sort_tag = barcode_tag,
                output_name = shard_prefix + ".consensus.bam"
        }

        if (emit_consensus_sorted) {
            call MergeFinalBams as MergeShardConsensusSortedBams {
                input:
                    bams = [select_first([LocalOverlapNonMito.consensus_sorted_bam]), select_first([LocalOverlapMito.consensus_sorted_bam])],
                    output_name = shard_prefix + ".consensus.sorted.bam"
            }
        }

        call CrossLocus {
            input:
                consensus_bam = MergeShardConsensusBams.output_bam,
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
            prune_pg_header = prune_pg_header_merge_final_bams,
            output_name = sample_name + ".merged.bam"
    }

    # Supplementary alignment BAMs are coordinate-sorted by construction for non-mito
    # and explicitly sorted in the mito branch before final merge.
    if (emit_supplementary_alignments) {
        Array[File] supplementary_bams_filtered = flatten([
            select_all(LocalOverlapNonMito.supplementary_alignments_bam),
            select_all(LocalOverlapMito.supplementary_alignments_bam)
        ])
        if (length(supplementary_bams_filtered) > 0) {
            call MergeFinalBams as MergeSupplementaryBams {
                input:
                    bams = supplementary_bams_filtered,
                    prune_pg_header = prune_pg_header_merge_final_bams,
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
        Array[File?] shard_consensus_sorted_bam = MergeShardConsensusSortedBams.merged_bam
        Array[File?] shard_consensus_sorted_bai = MergeShardConsensusSortedBams.merged_bai
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
        samtools index -@ 2 -o "~{basename(input_bam)}.bai" "~{input_bam}"
    >>>

    output {
        File bam_index = basename(input_bam) + ".bai"
    }

    runtime {
        cpu: 2
        memory: "2 GB"
        docker: docker_image
        disks: "local-disk ~{disk_gb} SSD"
        predefinedMachineType: "n2d-highcpu-2"
        preemptible: 3
    }
}

task LocalOverlap {
    input {
        File input_bam
        File input_bai
        File? reference_fasta
        File? reference_fai

        String output_prefix
        String barcode_tag
        String umi_tag
        String? tags
        Boolean no_consensus
        Boolean emit_supplementary_alignments
        Boolean emit_consensus_sorted
        Boolean mitochondrial_only = false
        String mitochondrial_contig_name = "chrM"

        Int? cpu
        Int? memory_gb

        String docker_image
    }

    Int task_cpu = select_first([cpu, 48])
    Int task_memory_gb = select_first([memory_gb, 48])
    Boolean use_predefined_machine_type = task_cpu == 2 || task_cpu == 4 || task_cpu == 8 || task_cpu == 16 || task_cpu == 32 || task_cpu == 48 || task_cpu == 64 || task_cpu == 80 || task_cpu == 96
    String machine_type = if (use_predefined_machine_type && task_memory_gb == task_cpu * 4)
        then "n2d-standard-${task_cpu}"
        else if (use_predefined_machine_type && task_memory_gb == task_cpu * 8)
        then "n2d-highmem-${task_cpu}"
        else if (use_predefined_machine_type && task_memory_gb == task_cpu)
        then "n2d-highcpu-${task_cpu}"
        else "n2d-custom-${task_cpu}-${task_memory_gb * 1024}"

    String tags_arg = if defined(tags) then "--tags " + tags else ""
    String supplementary_alignments_bam_path = output_prefix + ".supplementary_alignments.bam"
    String supplementary_alignments_arg = if emit_supplementary_alignments then "--sa-read-bam " + supplementary_alignments_bam_path else ""
    Int disk_gb = ceil(size(input_bam, "GB") * (if emit_consensus_sorted then 4 else 3)) + 20

    command <<<
        set -euo pipefail

        # Move input BAM and index to root working directory
        mv "~{input_bam}" "~{basename(input_bam)}"
        mv "~{input_bai}" "~{basename(input_bam)}.bai"

        chrom_list=$(samtools view -H "~{basename(input_bam)}" | awk -v mito="~{mitochondrial_contig_name}" -v mito_only="~{if mitochondrial_only then "1" else "0"}" '
            $1 == "@SQ" {
                contig = ""
                for (i = 2; i <= NF; ++i) {
                    if ($i ~ /^SN:/) {
                        contig = substr($i, 4)
                        break
                    }
                }
                if (contig == "") {
                    next
                }
                if ((mito_only == 1 && contig == mito) || (mito_only == 0 && contig != mito)) {
                    if (out != "") {
                        out = out " "
                    }
                    out = out contig
                }
            }
            END {
                print out
            }
        ')

        cudll_local_overlap \
            -i "~{basename(input_bam)}" \
            -o - \
            ~{if defined(reference_fasta) then "-r \"" + select_first([reference_fasta]) + "\"" else ""} \
            ~{tags_arg} \
            ~{supplementary_alignments_arg} \
            ~{if no_consensus then "--no-consensus" else ""} \
            -t ~{task_cpu} \
            -c "${chrom_list}" \
            --barcode ~{barcode_tag} \
            --umi ~{umi_tag} \
            --umi-hamming-only | \
            samtools sort --no-PG -@ ~{task_cpu} -t ~{barcode_tag} \
            -o "~{output_prefix}.consensus.bam" -

        if [ "~{emit_supplementary_alignments}" = "true" ] && [ "~{mitochondrial_only}" = "true" ]; then
            samtools sort --no-PG -@ ~{task_cpu} \
            -o "~{output_prefix}.supplementary_alignments.sorted.bam" \
            "~{supplementary_alignments_bam_path}"
            mv "~{output_prefix}.supplementary_alignments.sorted.bam" "~{supplementary_alignments_bam_path}"
        fi

        if [ "~{emit_consensus_sorted}" = "true" ]; then
            samtools sort --no-PG --write-index -@ ~{task_cpu} \
            -o "~{output_prefix}.consensus.sorted.bam##idx##~{output_prefix}.consensus.sorted.bam.bai" \
            "~{output_prefix}.consensus.bam"
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
        predefinedMachineType: "~{machine_type}"
        preemptible: 3
    }
}

task MergeTagSortedBams {
    input {
        Array[File] bams
        String sort_tag
        String output_name
    }

    Int diskGB = ceil(size(bams, "GB") + 20)

    command <<<
        set -euo pipefail

        samtools merge --no-PG -@ 2 -t "~{sort_tag}" -o "~{output_name}" ~{sep=' ' bams}
    >>>

    output {
        File output_bam = "~{output_name}"
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
    String machine_type = if defined(cpu) || defined(memory_gb) then "n2d-custom-${task_cpu}-${task_memory_gb * 1024}" else "n2d-highcpu-32"

    Int disk_gb = ceil(size(consensus_bam, "GB") * 3) + 20

    command <<<
        set -euo pipefail

        cudll_cross_locus \
            -i "~{consensus_bam}" \
            -o - \
            -t ~{task_cpu} \
            --barcode ~{barcode_tag} \
            --umi ~{umi_tag} \
            --identity ~{identity} \
            --umi-hamming-only \
            --rank-by-aligned-bases | \
            samtools sort --no-PG --write-index -@ ~{task_cpu} \
            -o "~{output_prefix}.consensus.homology_dedup.sorted.bam##idx##~{output_prefix}.consensus.homology_dedup.sorted.bam.bai" \
            -
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
        predefinedMachineType: "~{machine_type}"
        preemptible: 3
    }
}

task MergeFinalBams {
    input {
        Array[File] bams
        Boolean prune_pg_header = false
        String output_name
    }

    Int diskGB = ceil(size(bams, "GB") * (if prune_pg_header then 3.5 else 2.5) + 20)

    command <<<
        set -euo pipefail

        prune_pg_header() {
            local input_bam="$1"
            local output_bam="$2"
            local header_sam="$3"

            samtools view -H "${input_bam}" | awk '
                !/^@PG\t/ { print; next }
                /\tPN:minimap2(\t|$)/ { print; next }
                /\tPN:cudll_local_overlap(\t|$)/ {
                    if (local_line == "") {
                        local_line = $0
                        sub(/\tPP:[^\t]+/, "", local_line)
                    }
                    next
                }
                /\tPN:cudll_cross_locus(\t|$)/ {
                    if (cross_line == "") {
                        cross_line = $0
                        sub(/\tPP:[^\t]+/, "", cross_line)
                        sub(/\tPN:cudll_cross_locus/, "\tPN:cudll_cross_locus\tPP:cudll_local_overlap", cross_line)
                    }
                    next
                }
                { next }
                END {
                    if (local_line != "") print local_line
                    if (cross_line != "") print cross_line
                }
            ' > "${header_sam}"

            samtools reheader -P "${header_sam}" "${input_bam}" > "${output_bam}"
        }

        if [ "~{prune_pg_header}" = "true" ]; then
            declare -a merge_inputs=()

            for bam in ~{sep=' ' bams}; do
                pruned_bam="pruned_$(basename "$bam")"
                header_sam="${pruned_bam%.bam}.header.sam"

                prune_pg_header "$bam" "$pruned_bam" "$header_sam"
                merge_inputs+=("${pruned_bam}")
                rm -f "${header_sam}"
            done

            samtools merge --no-PG --write-index -p -@ 4 -o ~{output_name}##idx##~{output_name}.bai "${merge_inputs[@]}"
        else
            samtools merge --no-PG --write-index -p -@ 4 -o ~{output_name}##idx##~{output_name}.bai ~{sep=' ' bams}
        fi
    >>>

    output {
        File merged_bam = "~{output_name}"
        File merged_bai = "~{output_name}.bai"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        cpu: 4
        memory: "4 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 2
        predefinedMachineType: "n2d-highcpu-4"
    }
}
