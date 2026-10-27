version 1.0

# Pre- vs post-deduplication QC of aligned BAMs: collects the count tables (plots are made by Bam_Dedup_QC_Plot).
# Each module runs only if its inputs are given:
#   pre_bams (+ pre_bais)  -> scan each shard, merge  -> pre_* tables
#   post_bam (+ post_bai)  -> scan                    -> post_* tables
# Typical use: run once with the pre-dedup BAMs; later run with the dedup BAM only (the pre tables are not recomputed).
workflow Bam_Dedup_QC {
    input {
        String sample_name
        # Read end that position 0 refers to. "3p"/"5p" = library priming; use read_position_origin to override (e.g. priming "auto").
        String priming = "3p"
        String? read_position_origin

        Array[File]? pre_bams
        Array[File]? pre_bais
        File? post_bam
        File? post_bai

        Array[String] mito_contigs = ["chrM"]  # mitochondrial contig name(s) of the reference
        String cb_tag = "CB"
        String dc_tag = "DC"
        Int max_read_position = 10000
        Boolean include_supplementary = false
        String? exclude_tag  # skip reads that carry this BAM tag, e.g. "XR" (start-adapter rescued reads); needs a bam_dedup_qc image with --exclude-tag (0.2.0+)
        Int target_reads_per_chunk = 1000000
        Int scan_threads = 4
        String scan_machine_type = "c3d-highcpu-4"
        Int scan_extra_disk_gb = 10
        Int small_task_memory_gb = 8
        Int preemptible = 3
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/bam_dedup_qc:0.1.0"
    }

    Array[File] pre_bam_list = select_first([pre_bams, []])
    Array[File] pre_bai_list = select_first([pre_bais, []])
    Boolean run_pre = length(pre_bam_list) > 0
    Boolean run_post = defined(post_bam)

    # Fail before any BAM is localized (bad priming, missing/unpaired indexes, nothing to do).
    call CheckInputs {
        input:
            priming = priming, read_position_origin = read_position_origin,
            n_pre_bams = length(pre_bam_list), n_pre_bais = length(pre_bai_list),
            has_post_bam = defined(post_bam), has_post_bai = defined(post_bai),
            docker = docker, preemptible = preemptible
    }

    if (run_pre) {
        scatter (i in range(length(pre_bam_list))) {
            call Scan as PreScan {
                input:
                    bam = pre_bam_list[i], bai = pre_bai_list[i], label = "shard",
                    tar_only = true,
                    origin = CheckInputs.origin, mito_contigs = mito_contigs,
                    cb_tag = cb_tag,
                    dc_tag = dc_tag,
                    max_read_position = max_read_position,
                    include_supplementary = include_supplementary, exclude_tag = exclude_tag, target_reads_per_chunk = target_reads_per_chunk,
                    threads = scan_threads,
                    machine_type = scan_machine_type,
                    extra_disk_gb = scan_extra_disk_gb,
                    preemptible = preemptible,
                    docker = docker
            }
        }
        Int merge_chunk_size = 10
        Array[File] round0 = PreScan.tables_tar

        Int r1_total = length(round0)
        # round 1 only runs when the previous round left more than one chunk's worth of files
        if (r1_total > merge_chunk_size) {
            Int r1_num_chunks = (r1_total + merge_chunk_size - 1) / merge_chunk_size
            scatter (r1_idx in range(r1_num_chunks)) {
                Int r1_i0 = r1_idx * merge_chunk_size + 0
                Int r1_i1 = r1_idx * merge_chunk_size + 1
                Int r1_i2 = r1_idx * merge_chunk_size + 2
                Int r1_i3 = r1_idx * merge_chunk_size + 3
                Int r1_i4 = r1_idx * merge_chunk_size + 4
                Int r1_i5 = r1_idx * merge_chunk_size + 5
                Int r1_i6 = r1_idx * merge_chunk_size + 6
                Int r1_i7 = r1_idx * merge_chunk_size + 7
                Int r1_i8 = r1_idx * merge_chunk_size + 8
                Int r1_i9 = r1_idx * merge_chunk_size + 9
                Array[File] r1_chunk = flatten([
                    if (r1_i0 < r1_total) then [round0[r1_i0]] else [],
                    if (r1_i1 < r1_total) then [round0[r1_i1]] else [],
                    if (r1_i2 < r1_total) then [round0[r1_i2]] else [],
                    if (r1_i3 < r1_total) then [round0[r1_i3]] else [],
                    if (r1_i4 < r1_total) then [round0[r1_i4]] else [],
                    if (r1_i5 < r1_total) then [round0[r1_i5]] else [],
                    if (r1_i6 < r1_total) then [round0[r1_i6]] else [],
                    if (r1_i7 < r1_total) then [round0[r1_i7]] else [],
                    if (r1_i8 < r1_total) then [round0[r1_i8]] else [],
                    if (r1_i9 < r1_total) then [round0[r1_i9]] else []
                ])
                call MergeTarsPartial as r1_merge {
                    input:
                        label = sample_name + ".pre.r1." + r1_idx,
                        shard_tars = r1_chunk,
                        memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
                }
            }
        }
        Array[File] round1 = select_first([r1_merge.merged_tar, round0])

        Int r2_total = length(round1)
        # round 2 only runs when the previous round left more than one chunk's worth of files
        if (r2_total > merge_chunk_size) {
            Int r2_num_chunks = (r2_total + merge_chunk_size - 1) / merge_chunk_size
            scatter (r2_idx in range(r2_num_chunks)) {
                Int r2_i0 = r2_idx * merge_chunk_size + 0
                Int r2_i1 = r2_idx * merge_chunk_size + 1
                Int r2_i2 = r2_idx * merge_chunk_size + 2
                Int r2_i3 = r2_idx * merge_chunk_size + 3
                Int r2_i4 = r2_idx * merge_chunk_size + 4
                Int r2_i5 = r2_idx * merge_chunk_size + 5
                Int r2_i6 = r2_idx * merge_chunk_size + 6
                Int r2_i7 = r2_idx * merge_chunk_size + 7
                Int r2_i8 = r2_idx * merge_chunk_size + 8
                Int r2_i9 = r2_idx * merge_chunk_size + 9
                Array[File] r2_chunk = flatten([
                    if (r2_i0 < r2_total) then [round1[r2_i0]] else [],
                    if (r2_i1 < r2_total) then [round1[r2_i1]] else [],
                    if (r2_i2 < r2_total) then [round1[r2_i2]] else [],
                    if (r2_i3 < r2_total) then [round1[r2_i3]] else [],
                    if (r2_i4 < r2_total) then [round1[r2_i4]] else [],
                    if (r2_i5 < r2_total) then [round1[r2_i5]] else [],
                    if (r2_i6 < r2_total) then [round1[r2_i6]] else [],
                    if (r2_i7 < r2_total) then [round1[r2_i7]] else [],
                    if (r2_i8 < r2_total) then [round1[r2_i8]] else [],
                    if (r2_i9 < r2_total) then [round1[r2_i9]] else []
                ])
                call MergeTarsPartial as r2_merge {
                    input:
                        label = sample_name + ".pre.r2." + r2_idx,
                        shard_tars = r2_chunk,
                        memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
                }
            }
        }
        Array[File] round2 = select_first([r2_merge.merged_tar, round1])

        Int r3_total = length(round2)
        # round 3 only runs when the previous round left more than one chunk's worth of files
        if (r3_total > merge_chunk_size) {
            Int r3_num_chunks = (r3_total + merge_chunk_size - 1) / merge_chunk_size
            scatter (r3_idx in range(r3_num_chunks)) {
                Int r3_i0 = r3_idx * merge_chunk_size + 0
                Int r3_i1 = r3_idx * merge_chunk_size + 1
                Int r3_i2 = r3_idx * merge_chunk_size + 2
                Int r3_i3 = r3_idx * merge_chunk_size + 3
                Int r3_i4 = r3_idx * merge_chunk_size + 4
                Int r3_i5 = r3_idx * merge_chunk_size + 5
                Int r3_i6 = r3_idx * merge_chunk_size + 6
                Int r3_i7 = r3_idx * merge_chunk_size + 7
                Int r3_i8 = r3_idx * merge_chunk_size + 8
                Int r3_i9 = r3_idx * merge_chunk_size + 9
                Array[File] r3_chunk = flatten([
                    if (r3_i0 < r3_total) then [round2[r3_i0]] else [],
                    if (r3_i1 < r3_total) then [round2[r3_i1]] else [],
                    if (r3_i2 < r3_total) then [round2[r3_i2]] else [],
                    if (r3_i3 < r3_total) then [round2[r3_i3]] else [],
                    if (r3_i4 < r3_total) then [round2[r3_i4]] else [],
                    if (r3_i5 < r3_total) then [round2[r3_i5]] else [],
                    if (r3_i6 < r3_total) then [round2[r3_i6]] else [],
                    if (r3_i7 < r3_total) then [round2[r3_i7]] else [],
                    if (r3_i8 < r3_total) then [round2[r3_i8]] else [],
                    if (r3_i9 < r3_total) then [round2[r3_i9]] else []
                ])
                call MergeTarsPartial as r3_merge {
                    input:
                        label = sample_name + ".pre.r3." + r3_idx,
                        shard_tars = r3_chunk,
                        memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
                }
            }
        }
        Array[File] round3 = select_first([r3_merge.merged_tar, round2])

        Int r4_total = length(round3)
        # round 4 only runs when the previous round left more than one chunk's worth of files
        if (r4_total > merge_chunk_size) {
            Int r4_num_chunks = (r4_total + merge_chunk_size - 1) / merge_chunk_size
            scatter (r4_idx in range(r4_num_chunks)) {
                Int r4_i0 = r4_idx * merge_chunk_size + 0
                Int r4_i1 = r4_idx * merge_chunk_size + 1
                Int r4_i2 = r4_idx * merge_chunk_size + 2
                Int r4_i3 = r4_idx * merge_chunk_size + 3
                Int r4_i4 = r4_idx * merge_chunk_size + 4
                Int r4_i5 = r4_idx * merge_chunk_size + 5
                Int r4_i6 = r4_idx * merge_chunk_size + 6
                Int r4_i7 = r4_idx * merge_chunk_size + 7
                Int r4_i8 = r4_idx * merge_chunk_size + 8
                Int r4_i9 = r4_idx * merge_chunk_size + 9
                Array[File] r4_chunk = flatten([
                    if (r4_i0 < r4_total) then [round3[r4_i0]] else [],
                    if (r4_i1 < r4_total) then [round3[r4_i1]] else [],
                    if (r4_i2 < r4_total) then [round3[r4_i2]] else [],
                    if (r4_i3 < r4_total) then [round3[r4_i3]] else [],
                    if (r4_i4 < r4_total) then [round3[r4_i4]] else [],
                    if (r4_i5 < r4_total) then [round3[r4_i5]] else [],
                    if (r4_i6 < r4_total) then [round3[r4_i6]] else [],
                    if (r4_i7 < r4_total) then [round3[r4_i7]] else [],
                    if (r4_i8 < r4_total) then [round3[r4_i8]] else [],
                    if (r4_i9 < r4_total) then [round3[r4_i9]] else []
                ])
                call MergeTarsPartial as r4_merge {
                    input:
                        label = sample_name + ".pre.r4." + r4_idx,
                        shard_tars = r4_chunk,
                        memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
                }
            }
        }
        Array[File] round4 = select_first([r4_merge.merged_tar, round3])

        call Merge as PreMerge {
            input:
                label = sample_name + ".pre",
                shard_tars = round4,
                memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
        }
    }

    if (run_post) {
        call Scan as PostScan {
            input:
                bam = select_first([post_bam]),
                bai = select_first([post_bai]),
                label = sample_name + ".post",
                tar_only = false,
                origin = CheckInputs.origin,
                mito_contigs = mito_contigs,
                cb_tag = cb_tag,
                dc_tag = dc_tag,
                max_read_position = max_read_position,
                include_supplementary = include_supplementary, exclude_tag = exclude_tag,
                target_reads_per_chunk = target_reads_per_chunk,
                threads = scan_threads,
                machine_type = scan_machine_type,
                extra_disk_gb = scan_extra_disk_gb,
                preemptible = preemptible,
                docker = docker
        }
    }

    output {
        # Stage-prefixed so a pre-only run and a later post run never write the same table.
        # Same names as the Bam_Dedup_QC_Plot inputs.
        File? pre_aligned_length_tsv = PreMerge.aligned_length_tsv
        File? pre_dc_tsv = PreMerge.dc_tsv
        File? pre_error_by_position_tsv = PreMerge.error_by_position_tsv
        File? pre_read_error_rates_tsv = PreMerge.read_error_rates_tsv
        File? pre_cb_counts_tsv_gz = PreMerge.cb_counts_tsv_gz
        File? pre_contig_counts_tsv = PreMerge.contig_counts_tsv
        File? pre_summary_json = PreMerge.summary_json
        File? post_aligned_length_tsv = PostScan.aligned_length_tsv
        File? post_dc_tsv = PostScan.dc_tsv
        File? post_error_by_position_tsv = PostScan.error_by_position_tsv
        File? post_read_error_rates_tsv = PostScan.read_error_rates_tsv
        File? post_cb_counts_tsv_gz = PostScan.cb_counts_tsv_gz
        File? post_contig_counts_tsv = PostScan.contig_counts_tsv
        File? post_summary_json = PostScan.summary_json
    }
}

task CheckInputs {
    input {
        String priming
        String? read_position_origin
        Int n_pre_bams
        Int n_pre_bais
        Boolean has_post_bam
        Boolean has_post_bai
        String docker
        Int preemptible
    }
    String choice = select_first([read_position_origin, priming])

    command <<<
        set -euo pipefail
        err() { echo "ERROR: $1" >&2; exit 1; }
        case "~{choice}" in
            3p|5p) echo "~{choice}" > origin.txt ;;
            *) err "cannot derive the read position origin from '~{choice}' (priming='~{priming}'); set read_position_origin to 3p or 5p" ;;
        esac
        [ ~{n_pre_bams} -eq ~{n_pre_bais} ] || err "pre_bams has ~{n_pre_bams} files but pre_bais has ~{n_pre_bais}"
        [ "~{has_post_bam}" = "~{has_post_bai}" ] || err "post_bam and post_bai must be given together"
        if [ ~{n_pre_bams} -eq 0 ] && [ "~{has_post_bam}" = "false" ]; then
            err "nothing to do: give pre_bams+pre_bais and/or post_bam+post_bai"
        fi
    >>>

    output {
        String origin = read_string("origin.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        preemptible: preemptible
    }
}

# One pass over one BAM; writes the additive tables (see the tool README for the definitions).
# With tar_only the individual tables are not delivered (their outputs are null), only the tar.
task Scan {
    input {
        File bam
        File bai
        String label
        String origin
        Boolean tar_only  # deliver only the tar (for shards that are merged); otherwise also the individual tables
        Array[String] mito_contigs
        String cb_tag
        String dc_tag
        Int max_read_position
        Boolean include_supplementary
        String? exclude_tag
        Int target_reads_per_chunk
        Int threads
        String machine_type
        Int extra_disk_gb
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(size(bam, "GB") + size(bai, "GB") + extra_disk_gb)

    command <<<
        set -euo pipefail
        # The scanner finds the index next to the BAM: bring each pair into one directory under matching names.
        b=$(basename ~{bam}); i=$(basename ~{bai})
        if [ "$i" != "$b.bai" ] && [ "$i" != "${b%.bam}.bai" ]; then
            echo "ERROR: index $i does not belong to $b" >&2; exit 1
        fi
        mkdir -p in out
        ln -s ~{bam} "in/$b"
        ln -s ~{bai} "in/$b.bai"
        bamqc_scan --bam "in/$b" --label ~{label} --outdir out --threads ~{threads} --origin ~{origin} \
            --max-len ~{max_read_position} --cb-tag ~{cb_tag} --dc-tag ~{dc_tag} \
            --target-reads-per-chunk ~{target_reads_per_chunk} ~{sep=" " prefix("--mito-contig ", mito_contigs)} ~{if include_supplementary then "--include-supplementary" else ""} ~{if defined(exclude_tag) then "--exclude-tag " + select_first([exclude_tag]) else ""}
        tar -cf ~{label}.tables.tar -C out .
        ~{if tar_only then "rm -r out" else ""}
    >>>

    output {
        File? aligned_length_tsv = "out/~{label}.aligned_length.tsv"
        File? dc_tsv = "out/~{label}.dc.tsv"
        File? error_by_position_tsv = "out/~{label}.error_by_position.tsv"
        File? read_error_rates_tsv = "out/~{label}.read_error_rates.tsv"
        File? cb_counts_tsv_gz = "out/~{label}.cb_counts.tsv.gz"
        File? contig_counts_tsv = "out/~{label}.contig_counts.tsv"
        File? summary_json = "out/~{label}.summary.json"
        File tables_tar = "~{label}.tables.tar"  # the seven tables above in one file, so Merge localizes one object per shard
    }

    runtime {
        docker: docker
        predefinedMachineType: machine_type
        cpu: threads
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}

# Sums the tables of several scans (e.g. the shard BAMs of one library), given as one tar per scan.
task MergeTarsPartial {
    # Merges a chunk of shard (or already merged) table tars into one tar of the same layout, for the next round of the merge tree.
    input {
        String label
        Array[File] shard_tars
        Int memory_gb
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(size(shard_tars, "GB") * 3 + 10)

    command <<<
        set -euo pipefail
        mkdir -p shards out
        prefixes=()
        i=0
        for t in ~{sep=" " shard_tars}; do
            mkdir shards/s$i
            tar -xf "$t" -C shards/s$i
            sj=$(ls shards/s$i/*.summary.json)
            prefixes+=("${sj%.summary.json}")
            i=$((i+1))
        done
        bamqc_merge --label ~{label} --outdir out "${prefixes[@]}"
        tar -cf ~{label}.tables.tar -C out .
    >>>

    output {
        File merged_tar = "~{label}.tables.tar"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}

task Merge {
    input {
        String label
        Array[File] shard_tars
        Int memory_gb
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(size(shard_tars, "GB") * 3 + 10)

    command <<<
        set -euo pipefail
        mkdir -p shards out
        prefixes=()
        i=0
        for t in ~{sep=" " shard_tars}; do
            mkdir shards/s$i
            tar -xf "$t" -C shards/s$i
            # every shard was scanned with the same label; take the table prefix from its summary file
            sj=$(ls shards/s$i/*.summary.json)
            prefixes+=("${sj%.summary.json}")
            i=$((i+1))
        done
        bamqc_merge --label ~{label} --outdir out "${prefixes[@]}"
    >>>

    output {
        File aligned_length_tsv = "out/~{label}.aligned_length.tsv"
        File dc_tsv = "out/~{label}.dc.tsv"
        File error_by_position_tsv = "out/~{label}.error_by_position.tsv"
        File read_error_rates_tsv = "out/~{label}.read_error_rates.tsv"
        File cb_counts_tsv_gz = "out/~{label}.cb_counts.tsv.gz"
        File contig_counts_tsv = "out/~{label}.contig_counts.tsv"
        File summary_json = "out/~{label}.summary.json"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
