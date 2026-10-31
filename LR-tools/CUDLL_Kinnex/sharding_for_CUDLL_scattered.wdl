version 1.1

import "https://raw.githubusercontent.com/broadinstitute/MDL-workflows/main/misc_tools/extract_barcode_readcounts_from_aligned_bam/scBarcodeCounting.wdl" as counting
import "https://raw.githubusercontent.com/broadinstitute/MDL-workflows/main/SBX/shard_bams_by_barcode_group/shard_bams_by_barcode_group.wdl" as sharding
import "https://raw.githubusercontent.com/broadinstitute/MDL-workflows/main/SBX/CUDLL_scattered/CUDLL_scattered.wdl" as cudll

# Pipeline: (1) count reads per cell barcode across all BAMs
#           (2) detect knee, group barcodes into read-balanced groups, split/merge BAMs per group
#           (3) run CUDLL on the per-group BAMs as the "shards" of one sample
workflow sharding_for_CUDLL_scattered {
    input {
        Array[File]+ input_bams
        Array[File]+ input_bais
        String sample_name

        # shared
        String barcode_tag = "CB"

        # step 1: barcode counting
        String count_umi_tag = "XM"
        Boolean count_skip_umi = false
        Int count_cpu = 8
        String count_docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/brca/sc-barcode-counting:1.1"

        # step 2: sharding
        Int target_reads_per_group = 300000000
        Int min_rank_search = 2000
        Int max_rank_search = 200000
        Float slope_threshold = -2.0
        Float min_prominence = 0.8
        Boolean exclude_below_threshold = false

        # step 3: CUDLL
        String priming
        String docker_image_cudll
        String cudll_umi_tag = "UB"
        File? reference_fasta
        File? reference_fai
        Float identity = 0.95
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
        Int? local_overlap_mito_memory_gb
        Boolean local_overlap_mito_memory_auto = true
    }

    call counting.CountBarcodesWorkflow as CountBarcodes {
        input:
            bams = input_bams,
            bais = input_bais,
            cb_tag = barcode_tag,
            umi_tag = count_umi_tag,
            skip_umi = count_skip_umi,
            cpu = count_cpu,
            docker_image = count_docker_image
    }

    call sharding.Shard_Bams_By_Barcode_Group as ShardBams {
        input:
            input_bams = input_bams,
            input_bam_indexes = input_bais,
            barcode_tag = barcode_tag,
            counts_table = CountBarcodes.merged_counts,
            target_reads_per_group = target_reads_per_group,
            min_rank_search = min_rank_search,
            max_rank_search = max_rank_search,
            slope_threshold = slope_threshold,
            min_prominence = min_prominence,
            exclude_below_threshold = exclude_below_threshold
    }

    call cudll.CUDLL_scattered as CUDLL {
        input:
            input_bams = ShardBams.merged_group_bams,
            input_bais = ShardBams.merged_group_bam_indexes,
            mito_memory_gb_per_shard = ShardBams.group_mito_memory_gb,
            sample_name = sample_name,
            barcode_tag = barcode_tag,
            umi_tag = cudll_umi_tag,
            priming = priming,
            docker_image_cudll = docker_image_cudll,
            reference_fasta = reference_fasta,
            reference_fai = reference_fai,
            identity = identity,
            tags = tags,
            no_consensus = no_consensus,
            umi_allow_indel = umi_allow_indel,
            emit_supplementary_alignments = emit_supplementary_alignments,
            emit_consensus_sorted = emit_consensus_sorted,
            prune_pg_header_merge_final_bams = prune_pg_header_merge_final_bams,
            mitochondrial_contig_name = mitochondrial_contig_name,
            cpu = cpu,
            memory_gb = memory_gb,
            local_overlap_mito_cpu = local_overlap_mito_cpu,
            local_overlap_mito_memory_gb = local_overlap_mito_memory_gb,
            local_overlap_mito_memory_auto = local_overlap_mito_memory_auto
    }

    output {
        # step 1
        Array[File] per_bam_counts = CountBarcodes.counts
        File merged_counts = CountBarcodes.merged_counts
        # step 2
        File group_table = ShardBams.group_table
        Int n_groups = ShardBams.n_groups
        File group_memory_table = ShardBams.group_memory_table
        Array[File] merged_group_bams = ShardBams.merged_group_bams
        # step 3
        File CUDLL_bam = CUDLL.CUDLL_bam
        File CUDLL_bai = CUDLL.CUDLL_bai
        File? CUDLL_supplementary_bam = CUDLL.CUDLL_supplementary_bam
        File? CUDLL_supplementary_bai = CUDLL.CUDLL_supplementary_bai
        Array[File] CUDLL_pass2_shards_bam = CUDLL.CUDLL_pass2_shards_bam
        Array[File] CUDLL_pass2_shards_bai = CUDLL.CUDLL_pass2_shards_bai
        Array[File?] CUDLL_pass1_only_bam = CUDLL.CUDLL_pass1_only_bam
        Array[File?] CUDLL_pass1_only_bai = CUDLL.CUDLL_pass1_only_bai
        Int CUDLL_total_reads = CUDLL.CUDLL_total_reads
    }
}
