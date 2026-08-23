version 1.1

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
    }

    String excl_arg = if (exclude_below_threshold) then "--exclude-below-threshold" else ""
    Int diskGB = ceil(size(counts_table, "GB") * 2 + 10)

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
        cpu:         1
        memory:      "4 GB"
        disks:       "local-disk ~{diskGB} HDD"
        preemptible: 3
    }
}


task Merge_And_Split_Batch_Bams {
    input {
        Array[File] bams
        File        group_table
        Int         batch_index
        String      barcode_tag = "CB"
    }

    # Sized from 2026-08-22 sweep. c3d-highcpu-4 × batch=25 × --compression-level 1
    # gave $0.309/pipeline vs $0.498 for the previous n2d-highcpu-2 × batch=10 ×
    # default level=6 (38% cheaper compute). Efficiency: 87% at 4 threads on c3d.
    # At level=1 BGZF is 2-3x faster than the previous default level=6; outputs
    # are ~10-15% larger but go straight into Step 3's re-merge (also level=1).
    Int diskGB = ceil(size(bams, "GB") * 2 + 20)

    command <<<
        set -euo pipefail

        /usr/local/bin/split_bams_by_cb_group \
            --bams ~{sep=' ' bams} \
            --group-table ~{group_table} \
            --barcode-tag ~{barcode_tag} \
            --threads 4 \
            --compression-level 1 \
            --output-prefix batch~{batch_index}_group_
    >>>

    output {
        # glob returns files in alphabetical order: batch_group_0001.bam ... batch_group_NNNN.bam
        # stable ordering required for transpose() downstream.
        Array[File] group_bams = glob("batch~{batch_index}_group_*.bam")
    }

    runtime {
        docker:      "us-central1-docker.pkg.dev/methods-dev-lab/mdl-cudll/pysam-samtools:latest"
        cpu:         4
        memory:      "8 GB"
        predefinedMachineType: "c3d-highcpu-4"
        disks:       "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Merge_Group_Bams {
    input {
        Array[File] bams
        String      output_name = "merged.bam"
    }

    String output_index_name = "~{output_name}.bai"
    Int diskGB = ceil(size(bams, "GB") * 2.5 + 20)

    command <<<
        set -euo pipefail

        # k-way merge of coordinate-sorted inputs; output remains sorted.
        # level=1 output because these BAMs feed straight into CUDLL (which
        # decodes+re-encodes them anyway). Level 1 is ~3x faster than the
        # samtools default (level 6) with ~10-15% larger output.
        samtools merge --no-PG --write-index -p -@ 2 \
            --output-fmt-option level=1 \
            -o ~{output_name}##idx##~{output_index_name} \
            ~{sep=' ' bams}
    >>>

    output {
        File merged_bam = "~{output_name}"
        File merged_bam_index = "~{output_index_name}"
    }

    runtime {
        docker:      "us-central1-docker.pkg.dev/methods-dev-lab/samtools/samtools:latest"
        cpu:         2
        memory:      "2 GB"
        predefinedMachineType: "n2d-highcpu-2"
        disks:       "local-disk ~{diskGB} SSD"
        preemptible: 2
    }
}


workflow Shard_Bams_By_Barcode_Group {
    meta {
        description: "Detect the knee cutoff in a barcode counts table, assign barcodes to read-balanced groups, then split minimap2 BAM shards so that each output BAM contains all reads for one barcode group."
    }

    input {
        Array[File] input_bams
        Array[File] input_bam_indexes   # accepted for localization; not used directly

        String  barcode_tag = "CB"
        File    counts_table

        Int     target_reads_per_group = 300000000
        String  counts_column          = "post_count"
        String  barcode_column         = "barcode"
        Int     min_rank_search        = 2000
        Int     max_rank_search        = 200000
        Float   slope_threshold        = -2.0
        Float   min_prominence         = 0.8
        Boolean exclude_below_threshold = false
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
    }

    # ------------------------------------------------------------------
    # Step 2: Batch BAMs in groups of 25, then merge each batch in memory
    # while splitting directly by barcode group.
    #
    # batch_size=25 from the 2026-08-22 sweep: at c3d-highcpu-4 (8 GB RAM),
    # b=25 amortizes the ~12s barcode-map load and yields ~4-5% better $/GB
    # than b=10 while staying well under the RAM budget. Do NOT raise b above
    # 25 without moving to a larger-RAM machine — b=25 needs > 2 GB RAM.
    # ------------------------------------------------------------------
    Int batch_size   = 25
    Int total_bams   = length(input_bams)
    Int n_batches    = (total_bams + batch_size - 1) / batch_size

    scatter (b_idx in range(n_batches)) {
        # Inner scatter builds a File? per slot; outer collects into Array[File?]
        # which select_all trims. Cheaper than unrolling 25 conditional slots.
        scatter (j in range(batch_size)) {
            Int slot_idx = b_idx * batch_size + j
            File? maybe_bam = if (slot_idx < total_bams) then input_bams[slot_idx] else None
        }
        Array[File] batch_bams = select_all(maybe_bam)

        call Merge_And_Split_Batch_Bams {
            input:
                bams        = batch_bams,
                group_table = Assign_Barcode_Groups.group_table,
                batch_index = b_idx,
                barcode_tag = barcode_tag,
        }
    }

    # ------------------------------------------------------------------
    # Step 3: Transpose batch x group matrix -> group x batch, then merge
    # each group's BAMs into a single output BAM.
    #
    # transpose() requires every inner array to have the same length, which
    # is guaranteed because split_bams_by_cb_group always emits exactly
    # n_groups files (even if empty) with zero-padded names that glob in order.
    # ------------------------------------------------------------------
    Array[Array[File]] group_bam_matrix = transpose(Merge_And_Split_Batch_Bams.group_bams)

    scatter (g in zip(range(length(group_bam_matrix)), group_bam_matrix)) {
        call Merge_Group_Bams {
            input:
                bams        = g.right,
                output_name = "group_~{g.left}.bam",
        }
    }

    output {
        Array[File] merged_group_bams = Merge_Group_Bams.merged_bam  # group_0.bam, group_1.bam, ...
        Array[File] merged_group_bam_indexes = Merge_Group_Bams.merged_bam_index
        File        group_table       = Assign_Barcode_Groups.group_table
        Int         n_groups          = Assign_Barcode_Groups.n_groups
    }
}
