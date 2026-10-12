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
                    bam = pre_bam_list[i], bai = pre_bai_list[i], label = "shard", origin = CheckInputs.origin, mito_contigs = mito_contigs,
                    cb_tag = cb_tag, dc_tag = dc_tag, max_read_position = max_read_position,
                    include_supplementary = include_supplementary, target_reads_per_chunk = target_reads_per_chunk,
                    threads = scan_threads, machine_type = scan_machine_type, extra_disk_gb = scan_extra_disk_gb,
                    preemptible = preemptible, docker = docker
            }
        }
        call Merge as PreMerge {
            input:
                label = sample_name + ".pre",
                aligned_length_tsvs = PreScan.aligned_length_tsv, dc_tsvs = PreScan.dc_tsv,
                error_by_position_tsvs = PreScan.error_by_position_tsv, read_error_rates_tsvs = PreScan.read_error_rates_tsv, contig_counts_tsvs = PreScan.contig_counts_tsv, cb_counts_tsvs = PreScan.cb_counts_tsv_gz,
                summary_jsons = PreScan.summary_json,
                memory_gb = small_task_memory_gb, preemptible = preemptible, docker = docker
        }
    }

    if (run_post) {
        call Scan as PostScan {
            input:
                bam = select_first([post_bam]),
                bai = select_first([post_bai]),
                label = sample_name + ".post", origin = CheckInputs.origin, mito_contigs = mito_contigs,
                cb_tag = cb_tag, dc_tag = dc_tag, max_read_position = max_read_position,
                include_supplementary = include_supplementary, target_reads_per_chunk = target_reads_per_chunk,
                threads = scan_threads, machine_type = scan_machine_type, extra_disk_gb = scan_extra_disk_gb,
                preemptible = preemptible, docker = docker
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

# One pass over one BAM; writes the five additive tables (see the tool README for the definitions).
task Scan {
    input {
        File bam
        File bai
        String label
        String origin
        Array[String] mito_contigs
        String cb_tag
        String dc_tag
        Int max_read_position
        Boolean include_supplementary
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
            --target-reads-per-chunk ~{target_reads_per_chunk} ~{sep=" " prefix("--mito-contig ", mito_contigs)} ~{if include_supplementary then "--include-supplementary" else ""}
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
        predefinedMachineType: machine_type
        cpu: threads
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}

# Sums the tables of several scans (e.g. the shard BAMs of one library).
task Merge {
    input {
        String label
        Array[File] aligned_length_tsvs
        Array[File] dc_tsvs
        Array[File] error_by_position_tsvs
        Array[File] read_error_rates_tsvs
        Array[File] contig_counts_tsvs
        Array[File] cb_counts_tsvs
        Array[File] summary_jsons
        Int memory_gb
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(size(cb_counts_tsvs, "GB") * 4 + 10)

    command <<<
        set -euo pipefail
        mkdir -p shards out
        a=(~{sep=" " aligned_length_tsvs}); d=(~{sep=" " dc_tsvs}); e=(~{sep=" " error_by_position_tsvs}); r=(~{sep=" " read_error_rates_tsvs})
        c=(~{sep=" " cb_counts_tsvs}); k=(~{sep=" " contig_counts_tsvs}); s=(~{sep=" " summary_jsons})
        prefixes=()
        for i in "${!a[@]}"; do
            ln -s "${a[$i]}" shards/s$i.aligned_length.tsv
            ln -s "${d[$i]}" shards/s$i.dc.tsv
            ln -s "${e[$i]}" shards/s$i.error_by_position.tsv
            ln -s "${r[$i]}" shards/s$i.read_error_rates.tsv
            ln -s "${c[$i]}" shards/s$i.cb_counts.tsv.gz
            ln -s "${k[$i]}" shards/s$i.contig_counts.tsv
            ln -s "${s[$i]}" shards/s$i.summary.json
            prefixes+=("shards/s$i")
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
