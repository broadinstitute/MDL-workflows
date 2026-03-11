version 1.0

# Stub for Dockstore registration only.
# Replace with shard_bams_by_barcode_group.wdl.real to restore the real workflow.

task AssignBarcodeGroupsStub {
    input {
        File    counts_table
        String  counts_column          = "post_count"
        String  barcode_column         = "barcode"
        Int     min_rank_search        = 2000
        Int     max_rank_search        = 200000
        Float   slope_threshold        = -2.0
        Float   min_prominence         = 0.8
        Boolean exclude_below_threshold = false
        Int     target_reads_per_group = 300000000
    }
    command {}
    output {
        File group_table = counts_table
        Int  n_groups    = 1
    }
    runtime { docker: "ubuntu:24.04" }
}

workflow Shard_Bams_By_Barcode_Group {
    meta {
        description: "Detect the knee cutoff in a barcode counts table, assign barcodes to read-balanced groups, then split minimap2 BAM shards so that each output BAM contains all reads for one barcode group."
    }

    input {
        Array[File] input_bams
        Array[File] input_bam_indexes

        String  barcode_tag            = "CB"
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

    call AssignBarcodeGroupsStub {
        input:
            counts_table            = counts_table,
            counts_column           = counts_column,
            barcode_column          = barcode_column,
            min_rank_search         = min_rank_search,
            max_rank_search         = max_rank_search,
            slope_threshold         = slope_threshold,
            min_prominence          = min_prominence,
            exclude_below_threshold = exclude_below_threshold,
            target_reads_per_group  = target_reads_per_group
    }

    output {
        Array[File] merged_group_bams = input_bams
        File        group_table       = AssignBarcodeGroupsStub.group_table
        Int         n_groups          = AssignBarcodeGroupsStub.n_groups
    }
}
