version 1.0

workflow CUDLL_scattered {
    input {
        Array[File]+ input_bams
        Array[File]? input_bais
        File? reference_fasta
        File? reference_fai

        String sample_name
        String barcode_tag = "CB"
        String umi_tag = "UB"
        Float identity = 0.95
        String priming
        String? tags
        Boolean no_consensus = false
        Boolean umi_allow_indel = false
        Boolean emit_supplementary_alignments = true
        Boolean emit_consensus_sorted = false
        Boolean prune_pg_header_merge_final_bams = false
        String mitochondrial_contig_name = "chrM"

        Int? cpu
        Int? memory_gb
        Int? local_overlap_mito_cpu
        # Mito-pass memory, in order: local_overlap_mito_memory_gb (all shards) > auto (mito_memory_gb_per_shard, one per BAM) > old default.
        Int? local_overlap_mito_memory_gb
        Boolean local_overlap_mito_memory_auto = true
        Array[Int]? mito_memory_gb_per_shard

        String docker_image_cudll
    }

    Boolean mito_memory_auto_available = local_overlap_mito_memory_auto && defined(mito_memory_gb_per_shard) && length(select_first([mito_memory_gb_per_shard])) == length(input_bams)
    Boolean need_index = !defined(input_bais) || length(select_first([input_bais])) != length(input_bams)

    scatter (bam_idx in range(length(input_bams))) {
        File input_bam = input_bams[bam_idx]
        String shard_prefix = basename(input_bam, ".bam")

        if (need_index) {
            call CreateIndex { input: input_bam = input_bam, docker_image = docker_image_cudll }
        }

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
                priming = priming,
                tags = tags,
                no_consensus = no_consensus,
                umi_allow_indel = umi_allow_indel,
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
                priming = priming,
                tags = tags,
                no_consensus = no_consensus,
                umi_allow_indel = umi_allow_indel,
                emit_supplementary_alignments = emit_supplementary_alignments,
                emit_consensus_sorted = emit_consensus_sorted,
                mitochondrial_only = true,
                mitochondrial_contig_name = mitochondrial_contig_name,
                cpu = select_first([local_overlap_mito_cpu, 4]),
                memory_gb = if defined(local_overlap_mito_memory_gb) then local_overlap_mito_memory_gb else if mito_memory_auto_available then select_first([mito_memory_gb_per_shard])[bam_idx] else if defined(local_overlap_mito_cpu) then select_first([local_overlap_mito_cpu]) * 8 else if defined(memory_gb) then memory_gb * 8 else 32,
                docker_image = docker_image_cudll
        }

        call MergeTagSortedBams as MergeShardConsensusBams {
            input:
                bams = [LocalOverlapNonMito.consensus_bam, LocalOverlapMito.consensus_bam],
                sort_tag = barcode_tag,
                prune_pg_header = prune_pg_header_merge_final_bams,
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
                umi_allow_indel = umi_allow_indel,
                cpu = cpu,
                memory_gb = memory_gb,
                docker_image = docker_image_cudll
        }
    }

    call MergeFinalBams {
        input:
            bams = CrossLocus.final_bam,
            prune_pg_header = prune_pg_header_merge_final_bams,
            output_name = sample_name + ".CUDLL.bam"
    }

    # Supplementary BAMs: non-mito sorted by construction; mito branch sorts explicitly before merge.
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
        File   CUDLL_bam                       = MergeFinalBams.merged_bam
        File   CUDLL_bai                       = MergeFinalBams.merged_bai
        File?  CUDLL_supplementary_bam         = MergeSupplementaryBams.merged_bam
        File?  CUDLL_supplementary_bai         = MergeSupplementaryBams.merged_bai
        Array[File]  CUDLL_pass2_shards_bam    = CrossLocus.final_bam
        Array[File]  CUDLL_pass2_shards_bai    = CrossLocus.final_bai
        Array[File?] CUDLL_pass1_only_bam      = MergeShardConsensusSortedBams.merged_bam
        Array[File?] CUDLL_pass1_only_bai      = MergeShardConsensusSortedBams.merged_bai
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
        String priming
        String? tags
        Boolean no_consensus
        Boolean umi_allow_indel = false
        Boolean emit_supplementary_alignments
        Boolean emit_consensus_sorted
        Boolean mitochondrial_only = false
        String mitochondrial_contig_name = "chrM"

        Int? cpu
        Int? memory_gb

        String docker_image
    }

    Int task_cpu = select_first([cpu, 8])
    Int task_memory_gb = select_first([memory_gb, 32])
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int cpu_tier = if task_cpu <= 4 then 4
        else if task_cpu <= 8 then 8
        else if task_cpu <= 16 then 16
        else if task_cpu <= 30 then 30
        else if task_cpu <= 60 then 60
        else if task_cpu <= 90 then 90
        else if task_cpu <= 180 then 180
        else 360
    Int mem_tier = if task_memory_gb <= 32 then 4
        else if task_memory_gb <= 64 then 8
        else if task_memory_gb <= 128 then 16
        else if task_memory_gb <= 240 then 30
        else if task_memory_gb <= 480 then 60
        else if task_memory_gb <= 720 then 90
        else if task_memory_gb <= 1440 then 180
        else 360
    Int effective_cpu = if cpu_tier >= mem_tier then cpu_tier else mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int highcpu_ram = if effective_cpu == 4 then 8
        else if effective_cpu == 8 then 16
        else if effective_cpu == 16 then 32
        else if effective_cpu == 30 then 59
        else if effective_cpu == 60 then 118
        else if effective_cpu == 90 then 177
        else if effective_cpu == 180 then 354
        else 708
    String machine_type = if task_memory_gb <= highcpu_ram
        then "c3d-highcpu-${effective_cpu}"
        else if task_memory_gb <= effective_cpu * 4
        then "c3d-standard-${effective_cpu}"
        else "c3d-highmem-${effective_cpu}"

    String tags_arg = if defined(tags) then "--tags " + tags else ""
    String supplementary_alignments_bam_path = output_prefix + ".supplementary_alignments.bam"
    String supplementary_alignments_arg = if emit_supplementary_alignments then "--sa-read-bam " + supplementary_alignments_bam_path else ""
    String umi_flag = if umi_allow_indel then "--umi-allow-indel" else "--umi-hamming-only"
    Int disk_gb = ceil(size(input_bam, "GB") * (if emit_consensus_sorted then 4 else 3)) + 20

    command <<<
        set -euo pipefail

        mv "~{input_bam}" "~{basename(input_bam)}"
        mv "~{input_bai}" "~{basename(input_bam)}.bai"

        set +e
        _input_dev=$(df --output=source "$(dirname "$(readlink -f "~{basename(input_bam)}")")" 2>/dev/null | tail -n1)
        _base_dev=$(lsblk -no PKNAME "$_input_dev" 2>/dev/null | head -n1)
        if [ -n "$_base_dev" ] && [ -w "/sys/block/$_base_dev/queue/read_ahead_kb" ]; then
            echo 4096 > "/sys/block/$_base_dev/queue/read_ahead_kb"
            echo "[info] set read_ahead_kb=4096 on $_base_dev"
        else
            echo "[info] read_ahead tune skipped (base_dev='$_base_dev', input_dev='$_input_dev')"
        fi
        set -e

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

        if [ "~{mitochondrial_only}" = "true" ] && [ "~{priming}" != "auto" ]; then
            # --cb-sorted-input bounds memory by per-CB depth; skipped for priming=auto (calibration needs random access).
            samtools view -h -u -@ 1 "~{basename(input_bam)}" "${chrom_list}" | \
                samtools sort -u -@ 4 -t ~{barcode_tag} | \
                cudll_local_overlap \
                    -i - \
                    -o - \
                    ~{if defined(reference_fasta) then "-r \"" + select_first([reference_fasta]) + "\"" else ""} \
                    ~{tags_arg} \
                    ~{supplementary_alignments_arg} \
                    ~{if no_consensus then "--no-consensus" else ""} \
                    -t ~{effective_cpu} \
                    --barcode ~{barcode_tag} \
                    --umi ~{umi_tag} \
                    --priming ~{priming} \
                    ~{umi_flag} \
                    --cb-sorted-input | \
                samtools sort --no-PG -@ 4 -t ~{barcode_tag} \
                -o "~{output_prefix}.consensus.bam" -
        else
            cudll_local_overlap \
                -i "~{basename(input_bam)}" \
                -o - \
                ~{if defined(reference_fasta) then "-r \"" + select_first([reference_fasta]) + "\"" else ""} \
                ~{tags_arg} \
                ~{supplementary_alignments_arg} \
                ~{if no_consensus then "--no-consensus" else ""} \
                -t ~{effective_cpu} \
                -c "${chrom_list}" \
                --barcode ~{barcode_tag} \
                --umi ~{umi_tag} \
                --priming ~{priming} \
                ~{umi_flag} | \
                samtools sort --no-PG -@ 4 -t ~{barcode_tag} \
                -o "~{output_prefix}.consensus.bam" -
        fi

        if [ "~{emit_supplementary_alignments}" = "true" ] && [ "~{mitochondrial_only}" = "true" ]; then
            samtools sort --no-PG -@ 4 \
            -o "~{output_prefix}.supplementary_alignments.sorted.bam" \
            "~{supplementary_alignments_bam_path}"
            mv "~{output_prefix}.supplementary_alignments.sorted.bam" "~{supplementary_alignments_bam_path}"
        fi

        if [ "~{emit_consensus_sorted}" = "true" ]; then
            samtools sort --no-PG --write-index -@ 4 \
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

    # GCP Batch rejects predefinedMachineType + explicit cpu/memory together; omit them.
    runtime {
        docker: docker_image
        disks: "local-disk ~{disk_gb} SSD"
        predefinedMachineType: "~{machine_type}"
        preemptible: 3
        # Enables memoryRetryMultiplier: Cromwell only triggers a retry when maxRetries > 0.
        maxRetries: 2
    }
}

task MergeTagSortedBams {
    input {
        Array[File] bams
        String sort_tag
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

            samtools merge --no-PG -@ 2 -t "~{sort_tag}" -o "~{output_name}" "${merge_inputs[@]}"
        else
            samtools merge --no-PG -@ 2 -t "~{sort_tag}" -o "~{output_name}" ~{sep=' ' bams}
        fi
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
        Boolean umi_allow_indel = false

        Int? cpu
        Int? memory_gb

        String docker_image
    }

    Int task_cpu = select_first([cpu, 16])
    Int task_memory_gb = select_first([memory_gb, 16])
    # C3D has no custom shape; round up to the nearest fixed tier (4/8/16/30/60/90/180/360).
    Int cpu_tier = if task_cpu <= 4 then 4
        else if task_cpu <= 8 then 8
        else if task_cpu <= 16 then 16
        else if task_cpu <= 30 then 30
        else if task_cpu <= 60 then 60
        else if task_cpu <= 90 then 90
        else if task_cpu <= 180 then 180
        else 360
    Int mem_tier = if task_memory_gb <= 32 then 4
        else if task_memory_gb <= 64 then 8
        else if task_memory_gb <= 128 then 16
        else if task_memory_gb <= 240 then 30
        else if task_memory_gb <= 480 then 60
        else if task_memory_gb <= 720 then 90
        else if task_memory_gb <= 1440 then 180
        else 360
    Int effective_cpu = if cpu_tier >= mem_tier then cpu_tier else mem_tier
    # c3d-highcpu RAM is non-uniform per tier; use exact values.
    Int highcpu_ram = if effective_cpu == 4 then 8
        else if effective_cpu == 8 then 16
        else if effective_cpu == 16 then 32
        else if effective_cpu == 30 then 59
        else if effective_cpu == 60 then 118
        else if effective_cpu == 90 then 177
        else if effective_cpu == 180 then 354
        else 708
    String machine_type = if task_memory_gb <= highcpu_ram
        then "c3d-highcpu-${effective_cpu}"
        else if task_memory_gb <= effective_cpu * 4
        then "c3d-standard-${effective_cpu}"
        else "c3d-highmem-${effective_cpu}"

    String umi_flag = if umi_allow_indel then "--umi-allow-indel" else "--umi-hamming-only"
    Int disk_gb = ceil(size(consensus_bam, "GB") * 3) + 20

    command <<<
        set -euo pipefail

        cudll_cross_locus \
            -i "~{consensus_bam}" \
            -o - \
            -t ~{effective_cpu} \
            --barcode ~{barcode_tag} \
            --umi ~{umi_tag} \
            --identity ~{identity} \
            ~{umi_flag} \
            --rank-by-aligned-bases | \
            samtools sort --no-PG --write-index -@ 4 \
            -o "~{output_prefix}.consensus.homology_dedup.sorted.bam##idx##~{output_prefix}.consensus.homology_dedup.sorted.bam.bai" \
            -
    >>>

    output {
        File final_bam = "~{output_prefix}.consensus.homology_dedup.sorted.bam"
        File final_bai = "~{output_prefix}.consensus.homology_dedup.sorted.bam.bai"
    }

    # GCP Batch rejects predefinedMachineType + explicit cpu/memory together; omit them.
    runtime {
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
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 2
        predefinedMachineType: "c3d-highcpu-4"
    }
}
