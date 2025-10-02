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
        memory: "16 GB"
        cpu: 2
        disks: "local-disk 256 SSD"
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
        memory: "8 GB"
        cpu: 1
        disks: "local-disk 32 SSD"
    }
}

# Main workflow for processing a single sample
workflow add_metadata_to_ctat_fusion_call {
    input {
        File bam_file
        File bam_index
        File fusion_tsv
        String pool_name
        File metadata_tsv
        File path_extract_CB_py
        File path_add_metadata_py
        String docker_image = "python:3.9-slim"
    }

    call extract_barcodes_and_update_fusion {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            fusion_tsv = fusion_tsv,
            path_extract_CB_py = path_extract_CB_py,
            pool_name = pool_name,
            docker_image = docker_image
    }

    call add_metadata_to_fusion {
        input:
            fusion_with_tags_tsv = extract_barcodes_and_update_fusion.updated_fusion_tsv,
            metadata_tsv = metadata_tsv,
            path_add_metadata_py = path_add_metadata_py,
            pool_name = pool_name,
            docker_image = docker_image
    }

    output {
        File fusion_with_tags_tsv = extract_barcodes_and_update_fusion.updated_fusion_tsv
        File final_fusion_tsv = add_metadata_to_fusion.fusion_with_metadata
    }
}
