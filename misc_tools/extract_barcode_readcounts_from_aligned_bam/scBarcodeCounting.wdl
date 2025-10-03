version 1.0

workflow BarcodeCountingWorkflow {
    input {
        Array[File] bam_files
        Array[File] bam_indices  # .bai files
        File metadata_tsv
        File process_barcodes_script  # The process_barcodes.py script
        File merge_results_script  # The merge_results.py script
        String barcode_tag = "CB"
        String docker_image = "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-python:1.0.15"
    }

    # Process each BAM file in parallel
    scatter (i in range(length(bam_files))) {
        call ProcessSingleBAM {
            input:
                bam_file = bam_files[i],
                bam_index = bam_indices[i],
                metadata_tsv = metadata_tsv,
                process_barcodes_script = process_barcodes_script,
                barcode_tag = barcode_tag,
                docker_image = docker_image
        }
    }

    # Merge all results
    call MergeResults {
        input:
            patient_stats_files = ProcessSingleBAM.patient_stats,
            pool_barcode_files = ProcessSingleBAM.pool_barcode_details,
            unmatched_files = ProcessSingleBAM.unmatched_barcodes,
            merge_results_script = merge_results_script,
            docker_image = docker_image
    }

    output {
        File final_patient_summary = MergeResults.patient_summary
        File final_unmatched = MergeResults.unmatched_summary
        File final_all_barcodes = MergeResults.all_barcodes_combined
        Array[File] individual_patient_stats = ProcessSingleBAM.patient_stats
        Array[File] pool_barcode_details = ProcessSingleBAM.pool_barcode_details
        Array[File] individual_unmatched = ProcessSingleBAM.unmatched_barcodes
    }
}

task ProcessSingleBAM {
    input {
        File bam_file
        File bam_index
        File metadata_tsv
        File process_barcodes_script
        String barcode_tag
        String docker_image
        Int disk_size = 150
        Int memory_gb = 16
        Int cpu = 4
    }

    String bam_basename = basename(bam_file, ".bam")

    command <<<
        set -euo pipefail

        # Install dependencies
        pip install pysam pandas --quiet

        # Make script executable
        chmod +x ~{process_barcodes_script}

        # Run the Python script
        python ~{process_barcodes_script} \
            ~{bam_file} \
            ~{metadata_tsv} \
            --barcode-tag ~{barcode_tag} \
            --output-prefix ~{bam_basename}
    >>>

    output {
        File patient_stats = "~{bam_basename}_patient_stats.tsv"
        File pool_barcode_details = "~{bam_basename}_barcode_details.tsv"
        File unmatched_barcodes = "~{bam_basename}_unmatched.tsv"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        cpu: cpu
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 2
        maxRetries: 1
    }
}

task MergeResults {
    input {
        Array[File] patient_stats_files
        Array[File] pool_barcode_files
        Array[File] unmatched_files
        File merge_results_script
        String docker_image
        Int disk_size = 50
        Int memory_gb = 8
    }

    command <<<
        set -euo pipefail

        pip install pandas --quiet

        # Make script executable
        chmod +x ~{merge_results_script}

        # Run the merge script
        python ~{merge_results_script} \
            --patient-stats ~{sep=' ' patient_stats_files} \
            --barcode-details ~{sep=' ' pool_barcode_files} \
            --unmatched ~{sep=' ' unmatched_files} \
            --output-prefix merged
    >>>

    output {
        File patient_summary = "merged_patient_summary.tsv"
        File unmatched_summary = "merged_unmatched.tsv"
        File all_barcodes_combined = "merged_all_barcodes.tsv"
    }

    runtime {
        docker: docker_image
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        preemptible: 2
    }
}
