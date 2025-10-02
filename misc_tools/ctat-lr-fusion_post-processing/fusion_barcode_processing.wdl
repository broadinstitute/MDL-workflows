version 1.0

# Task 1: Extract cell barcodes and UMI from BAM files and update fusion TSV
task extract_barcodes_and_update_fusion {
    input {
        File bam_file
        File bam_index
        File fusion_tsv
        File path_extract_CB_py
        String pool_name
        String docker_image = "python:3.9-slim"
    }

    String output_filename = "${pool_name}.fusion_calls.with_tags.tsv"

    command <<<
        # Install required Python packages
        pip install pysam pandas

        # Copy the script to working directory
        cp ~{path_extract_CB_py} extract_tags.py

        # Run the Python script with arguments
        python extract_tags.py "~{bam_file}" "~{fusion_tsv}" "~{output_filename}"
    >>>

    output {
        File updated_fusion_tsv = output_filename
    }

    runtime {
        docker: docker_image
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 50 SSD"
    }
}

# Task 2: Add metadata fields to fusion calls based on cell barcodes
task add_metadata_to_fusion {
    input {
        File fusion_with_tags_tsv
        File metadata_tsv
        File path_add_metadata_py
        String pool_name
        String docker_image = "python:3.9-slim"
    }

    String output_filename = "${pool_name}.fusion_calls.with_tags.meta.tsv"

    command <<<
        # Install required Python packages
        pip install pandas

        # Copy the script to working directory
        cp ~{path_add_metadata_py} add_metadata.py

        # Run the Python script with arguments
        python add_metadata.py "~{fusion_with_tags_tsv}" "~{metadata_tsv}" "~{output_filename}"
    >>>

    output {
        File fusion_with_metadata = output_filename
    }

    runtime {
        docker: docker_image
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 20 SSD"
    }
}

# Main workflow that combines both tasks
workflow fusion_barcode_processing_workflow {
    input {
        Array[File] bam_files
        Array[File] bam_indices
        Array[File] fusion_tsvs
        Array[String] pool_names
        File metadata_tsv
        File path_extract_CB_py
        File path_add_metadata_py
        String docker_image = "python:3.9-slim"
    }

    # Process each pool
    scatter (i in range(length(bam_files))) {
        call extract_barcodes_and_update_fusion {
            input:
                bam_file = bam_files[i],
                bam_index = bam_indices[i],
                fusion_tsv = fusion_tsvs[i],
                path_extract_CB_py = path_extract_CB_py,
                pool_name = pool_names[i],
                docker_image = docker_image
        }

        call add_metadata_to_fusion {
            input:
                fusion_with_tags_tsv = extract_barcodes_and_update_fusion.updated_fusion_tsv,
                metadata_tsv = metadata_tsv,
                path_add_metadata_py = path_add_metadata_py,
                pool_name = pool_names[i],
                docker_image = docker_image
        }
    }

    output {
        Array[File] fusion_with_tags_tsvs = extract_barcodes_and_update_fusion.updated_fusion_tsv
        Array[File] final_fusion_tsvs = add_metadata_to_fusion.fusion_with_metadata
    }
}
