version 1.0

workflow LRAA_PostProcessing {
    input {
        String sample_name
        File LRAA_tracking_file      # Output from LRAA.wdl: mergedQuantTracking
        File LRAA_gtf_file          # Output from LRAA.wdl: mergedGTF
        File reference_gtf_file     # GENCODE reference annotation GTF
        
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
        Int memoryGB = 128
        Int diskSizeGB = 1024
        Int numThreads = 16
    }

    # Step 1: Convert single cell tracking to sparse matrix
    call singlecell_tracking_to_sparse_matrix {
        input:
            tracking_file = LRAA_tracking_file,
            output_prefix = sample_name,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    # Step 2: Use gffcompare to map LRAA isoform structures to reference annotation
    call gffcompare_mapping {
        input:
            reference_gtf = reference_gtf_file,
            lraa_gtf = LRAA_gtf_file,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    # Step 3: Update LRAA gene symbols using reference annotation
    call incorporate_gene_symbols {
        input:
            lraa_gtf = LRAA_gtf_file,
            reference_gtf = reference_gtf_file,
            gffcompare_tracking = gffcompare_mapping.tracking_file,
            id_mappings_file = singlecell_tracking_to_sparse_matrix.id_mappings_file,
            gene_sparseM_dir = singlecell_tracking_to_sparse_matrix.gene_sparseM_dir,
            isoform_sparseM_dir = singlecell_tracking_to_sparse_matrix.isoform_sparseM_dir,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    output {
        # Sparse matrix outputs from step 1
        File gene_cell_counts = singlecell_tracking_to_sparse_matrix.gene_cell_counts
        File isoform_cell_counts = singlecell_tracking_to_sparse_matrix.isoform_cell_counts
        File splice_pattern_cell_counts = singlecell_tracking_to_sparse_matrix.splice_pattern_cell_counts
        File id_mappings = singlecell_tracking_to_sparse_matrix.id_mappings_file
        File gene_sparseM_dir = singlecell_tracking_to_sparse_matrix.gene_sparseM_dir
        File isoform_sparseM_dir = singlecell_tracking_to_sparse_matrix.isoform_sparseM_dir
        File splice_pattern_sparseM_dir = singlecell_tracking_to_sparse_matrix.splice_pattern_sparseM_dir
        
        # gffcompare outputs from step 2
        File gffcompare_annotated_gtf = gffcompare_mapping.annotated_gtf
        File gffcompare_tracking = gffcompare_mapping.tracking_file
        File gffcompare_loci = gffcompare_mapping.loci_file
        File gffcompare_stats = gffcompare_mapping.stats_file
        
        # Annotated outputs from step 3
        File updated_gene_sparseM_dir = incorporate_gene_symbols.updated_gene_sparseM_dir
        File updated_isoform_sparseM_dir = incorporate_gene_symbols.updated_isoform_sparseM_dir
        File annotated_id_mappings = incorporate_gene_symbols.annotated_id_mappings
        File? updated_lraa_gtf = incorporate_gene_symbols.updated_lraa_gtf
    }
}

task singlecell_tracking_to_sparse_matrix {
    input {
        File tracking_file
        String output_prefix
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -euo pipefail
        
        # Decompress tracking file if it's gzipped
        if [[ "~{tracking_file}" == *.gz ]]; then
            gunzip -c ~{tracking_file} > tracking_decompressed.tsv
            tracking_input="tracking_decompressed.tsv"
        else
            tracking_input="~{tracking_file}"
        fi
        
        # Run the R script
        /usr/local/src/LRAA/util/sc/singlecell_tracking_to_sparse_matrix.R \
            --tracking "$tracking_input" \
            --output_prefix ~{output_prefix}
        
        # Create tar archives for the sparse matrix directories
        tar -czf ~{output_prefix}^gene-sparseM.tar.gz ~{output_prefix}^gene-sparseM/
        tar -czf ~{output_prefix}^isoform-sparseM.tar.gz ~{output_prefix}^isoform-sparseM/
        tar -czf ~{output_prefix}^splice_pattern-sparseM.tar.gz ~{output_prefix}^splice_pattern-sparseM/
    >>>

    output {
        File gene_cell_counts = "~{output_prefix}.gene_cell_counts.tsv"
        File isoform_cell_counts = "~{output_prefix}.isoform_cell_counts.tsv"
        File splice_pattern_cell_counts = "~{output_prefix}.splice_pattern_cell_counts.tsv"
        File id_mappings_file = "~{output_prefix}.gene_transcript_splicehashcode.tsv"
        File gene_sparseM_dir = "~{output_prefix}^gene-sparseM.tar.gz"
        File isoform_sparseM_dir = "~{output_prefix}^isoform-sparseM.tar.gz"
        File splice_pattern_sparseM_dir = "~{output_prefix}^splice_pattern-sparseM.tar.gz"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task gffcompare_mapping {
    input {
        File reference_gtf
        File lraa_gtf
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -euo pipefail
        
        # Run gffcompare
        gffcompare -r ~{reference_gtf} ~{lraa_gtf}
    >>>

    output {
        File annotated_gtf = "gffcmp.annotated.gtf"
        File tracking_file = "gffcmp.tracking"
        File loci_file = "gffcmp.loci"
        File stats_file = "gffcmp.stats"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}

task incorporate_gene_symbols {
    input {
        File lraa_gtf
        File reference_gtf
        File gffcompare_tracking
        File id_mappings_file
        File gene_sparseM_dir
        File isoform_sparseM_dir
        String docker
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        set -euo pipefail
        
        # Extract sparse matrix directories
        tar -xzf ~{gene_sparseM_dir}
        tar -xzf ~{isoform_sparseM_dir}
        
        # Get the directory names (remove .tar.gz extension and path)
        gene_dir=$(basename ~{gene_sparseM_dir} .tar.gz)
        isoform_dir=$(basename ~{isoform_sparseM_dir} .tar.gz)
        
        # Run the Python script to incorporate gene symbols
        /usr/local/src/LRAA/util/sc/incorporate_gene_symbols_in_sc_features.py \
            --LRAA_gtf ~{lraa_gtf} \
            --ref_gtf ~{reference_gtf} \
            --gffcompare_tracking ~{gffcompare_tracking} \
            --id_mappings ~{id_mappings_file} \
            --sparseM_dirs "$gene_dir" "$isoform_dir"
        
        # Create tar archives for the updated sparse matrix directories
        tar -czf updated_${gene_dir}.tar.gz "$gene_dir"/
        tar -czf updated_${isoform_dir}.tar.gz "$isoform_dir"/
    >>>

    output {
        File updated_gene_sparseM_dir = glob("updated_*gene-sparseM.tar.gz")[0]
        File updated_isoform_sparseM_dir = glob("updated_*isoform-sparseM.tar.gz")[0]
        File annotated_id_mappings = "~{id_mappings_file}.wAnnotIDs"
        File? updated_lraa_gtf = "~{lraa_gtf}.updated.gtf"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }
}
