version 1.1

# ---------------------------------------------------------------------------
# Shard BAMs by barcode group
#
# 1. Runs assign_barcode_groups.py on a whitelist counts table to detect the
#    first knee and assign each barcode to one of N groups (~300 M reads each).
#
# 2. Scatters over batches of 10 input BAMs.  Each batch is split into N
#    per-group BAMs (one per CB group) using split_bams_by_cb_group.py,
#    then each slice is coordinate-sorted in place.
#
# 3. Transposes the batch x group matrix and scatters over groups: all batch
#    outputs for the same group are merged into a single final sorted BAM.
#
# Motivation: 1k+ minimap2 shards shouldn't be merged+sorted as a single job 
# (very long running, limited to 1 core for the merge) and downstream deduplication 
# will have too much depth in highly expressed regions (CPU-bound, memory-intensive).
# Batching by 10 reduces task count to ~100+ and lets the final per-group merges
# run in parallel across N jobs.
# ---------------------------------------------------------------------------

task Assign_Barcode_Groups {
    input {
        File counts_table

        String counts_column         = "post_count"
        String barcode_column        = "barcode"
        Int    min_rank_search       = 2000
        Int    max_rank_search       = 200000
        Float  slope_threshold       = -2.0
        Float  min_prominence        = 0.8
        Boolean exclude_below_threshold = false
        Int    target_reads_per_group = 300000000

        Int    cpu      = 1
        Int    mem_gb   = 4
        Int?   disk_size_gb
    }

    String excl_arg = if exclude_below_threshold then "--exclude-below-threshold" else ""
    Int diskGB = select_first([disk_size_gb, ceil(size(counts_table, "GB") * 2 + 10)])

    command <<<
        set -euo pipefail

        python3 /usr/local/bin/assign_barcode_groups.py ~{counts_table} \
            --counts-column ~{counts_column} \
            --barcode-column ~{barcode_column} \
            --min-rank-search ~{min_rank_search} \
            --max-rank-search ~{max_rank_search} \
            --slope-threshold ~{slope_threshold} \
            --min-prominence ~{min_prominence} \
            ~{excl_arg} \
            --target-reads-per-group ~{target_reads_per_group} \
            --group-output group_table.tsv

        # Extract N (max group number) from the output table.
        awk -F'\t' 'NR > 1 && $2 + 0 > max { max = $2 + 0 } END { print max }' \
            group_table.tsv > n_groups.txt
    >>>

    output {
        File group_table = "group_table.tsv"
        Int  n_groups    = read_int("n_groups.txt")
    }

    runtime {
        docker:      "us-central1-docker.pkg.dev/methods-dev-lab/mdl-cudll/assign-barcode-groups:latest"
        cpu:         cpu
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{diskGB} HDD"
        preemptible: 3
    }
}


task Merge_And_Split_Batch {
    input {
        Array[File] bams
        File        group_table
        Int         n_groups
        Int         batch_index
        String      barcode_tag = "CB"

        Int    cpu       = 2
        Int    mem_gb    = 4
        Int?   disk_size_gb
    }

    Int diskGB = select_first([disk_size_gb, ceil(size(bams, "GB") * 2.5 + 20)])

    command <<<
        set -euo pipefail

        python3 /usr/local/bin/split_bams_by_cb_group.py \
            --bams ~{sep=' ' bams} \
            --group-table ~{group_table} \
            --n-groups ~{n_groups} \
            --barcode-tag ~{barcode_tag} \
            --output-prefix batch~{batch_index}_group_

        # Coordinate-sort each group BAM in place.  Run cpu single-threaded
        # sorts in parallel - more efficient than one multi-threaded sort for
        # small per-batch slices.
        printf '%s\n' batch~{batch_index}_group_*.bam | \
            xargs -P ~{cpu} -I{} bash -c 'samtools sort -o "${1%.bam}.sorted.bam" "$1"' _ {}
    >>>

    output {
        # glob returns files in sorted order; batch_group_0001.sorted.bam ... batch_group_NNNN.sorted.bam
        # gives the stable ordering required for transpose() downstream.
        Array[File] group_bams = glob("batch~{batch_index}_group_*.sorted.bam")
    }

    runtime {
        docker:      "us-central1-docker.pkg.dev/methods-dev-lab/mdl-cudll/pysam-samtools:latest"
        cpu:         cpu
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Merge_Group_Bams {
    input {
        Array[File] bams
        String      output_name = "merged.bam"

        Int    cpu     = 2
        Int    mem_gb  = 4
        Int?   disk_size_gb
    }

    Int diskGB = select_first([disk_size_gb, ceil(size(bams, "GB") * 2.5 + 20)])

    command <<<
        set -euo pipefail

        # k-way merge of coordinate-sorted inputs; output remains sorted.
        samtools merge -@ ~{cpu} -o ~{output_name} ~{sep=' ' bams}
    >>>

    output {
        File merged_bam = "~{output_name}"
    }

    runtime {
        docker:      "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        cpu:         cpu
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{diskGB} SSD"
        preemptible: 2
    }
}


workflow Shard_Bams_By_Barcode_Group {
    meta {
        description: "Detect the knee cutoff in a barcode counts table, assign barcodes to read-balanced groups, then split minimap2 BAM shards so that each output BAM contains all reads for one barcode group."
    }

    input {
        # ---- knee detection & group assignment ----
        File    counts_table
        String  counts_column          = "post_count"
        String  barcode_column         = "barcode"
        Int     min_rank_search        = 2000
        Int     max_rank_search        = 200000
        Float   slope_threshold        = -2.0
        Float   min_prominence         = 0.8
        Boolean exclude_below_threshold = false
        Int     target_reads_per_group = 300000000

        # ---- BAM inputs ----
        Array[File] input_bams
        Array[File] input_bam_indexes   # accepted as input for localization; not used directly

        # ---- barcode tag ----
        String  barcode_tag = "CB"

        # ---- resource overrides ----
        Int    assign_groups_cpu            = 1
        Int    assign_groups_mem_gb         = 4
        Int    split_batch_cpu              = 2
        Int    split_batch_mem_gb           = 4
        Int?   split_batch_disk_size_gb
        Int    merge_group_cpu              = 2
        Int    merge_group_mem_gb           = 4
        Int?   merge_group_disk_size_gb
    }

    # ------------------------------------------------------------------
    # Step 1: Detect first knee and assign each barcode to a read group.
    # ------------------------------------------------------------------
    call Assign_Barcode_Groups {
        input:
            counts_table             = counts_table,
            counts_column            = counts_column,
            barcode_column           = barcode_column,
            min_rank_search          = min_rank_search,
            max_rank_search          = max_rank_search,
            slope_threshold          = slope_threshold,
            min_prominence           = min_prominence,
            exclude_below_threshold  = exclude_below_threshold,
            target_reads_per_group   = target_reads_per_group,
            cpu                      = assign_groups_cpu,
            mem_gb                   = assign_groups_mem_gb
    }

    # ------------------------------------------------------------------
    # Step 2: Batch BAMs in groups of 10, then split each batch by group.
    #
    # Uses the select_all / None pattern to handle the last partial batch.
    # ------------------------------------------------------------------
    Int batch_size   = 10
    Int total_bams   = length(input_bams)
    Int n_batches    = (total_bams + batch_size - 1) / batch_size

    scatter (b_idx in range(n_batches)) {
        Int b0 = b_idx * batch_size + 0
        Int b1 = b_idx * batch_size + 1
        Int b2 = b_idx * batch_size + 2
        Int b3 = b_idx * batch_size + 3
        Int b4 = b_idx * batch_size + 4
        Int b5 = b_idx * batch_size + 5
        Int b6 = b_idx * batch_size + 6
        Int b7 = b_idx * batch_size + 7
        Int b8 = b_idx * batch_size + 8
        Int b9 = b_idx * batch_size + 9

        Array[File] batch_bams = select_all([
            if (b0 < total_bams) then input_bams[b0] else None,
            if (b1 < total_bams) then input_bams[b1] else None,
            if (b2 < total_bams) then input_bams[b2] else None,
            if (b3 < total_bams) then input_bams[b3] else None,
            if (b4 < total_bams) then input_bams[b4] else None,
            if (b5 < total_bams) then input_bams[b5] else None,
            if (b6 < total_bams) then input_bams[b6] else None,
            if (b7 < total_bams) then input_bams[b7] else None,
            if (b8 < total_bams) then input_bams[b8] else None,
            if (b9 < total_bams) then input_bams[b9] else None
        ])

        call Merge_And_Split_Batch {
            input:
                bams        = batch_bams,
                group_table = Assign_Barcode_Groups.group_table,
                n_groups    = Assign_Barcode_Groups.n_groups,
                batch_index = b_idx,
                barcode_tag = barcode_tag,
                cpu         = split_batch_cpu,
                mem_gb      = split_batch_mem_gb,
                disk_size_gb = split_batch_disk_size_gb
        }
    }

    # ------------------------------------------------------------------
    # Step 3: Transpose batch x group matrix -> group x batch, then merge
    # each group's BAMs into a single output BAM.
    #
    # transpose() requires every inner array to have the same length, which
    # is guaranteed because split_bams_by_cb_group.py always emits exactly
    # n_groups files (even if empty) with zero-padded names that glob in order.
    # ------------------------------------------------------------------
    Array[Array[File]] per_group_bams = transpose(Merge_And_Split_Batch.group_bams)

    scatter (group_bams in per_group_bams) {
        call Merge_Group_Bams {
            input:
                bams         = group_bams,
                cpu          = merge_group_cpu,
                mem_gb       = merge_group_mem_gb,
                disk_size_gb = merge_group_disk_size_gb
        }
    }

    output {
        Array[File] merged_group_bams = Merge_Group_Bams.merged_bam
        File        group_table       = Assign_Barcode_Groups.group_table
        Int         n_groups          = Assign_Barcode_Groups.n_groups
    }
}
