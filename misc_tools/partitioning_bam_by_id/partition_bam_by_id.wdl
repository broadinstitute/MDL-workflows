version 1.0

## Workflow to partition BAM files by cell clusters

workflow partition_bam_by_cell_cluster_workflow {
    input {
        String sample_id
        String partitioning_id
        File inputBAM
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }
    
    # Main partitioning task
    call partition_bam_by_cell_cluster {
        input:
            sample_id = sample_id,
            partitioning_id = partitioning_id,
            inputBAM = inputBAM,
            docker = docker
    }
    
    output {
        Array[File] partitioned_bams = partition_bam_by_cell_cluster.partitioned_bams
    }
    
    meta {
        description: "Partition BAM file by cell clusters"
        author: "Your Name"
        version: "1.0"
    }
    
    parameter_meta {
        sample_id: "Unique identifier for the sample"
        partitioning_id: "Partitioning identifier for cell cluster assignments"
        inputBAM: "Input BAM file to partition"
        docker: "Docker image containing partition_bam_by_cell_cluster.py"
        partitioned_bams: "Array of partitioned BAM files, one per cluster"
    }
}

# Main task to partition BAM by cell clusters
task partition_bam_by_cell_cluster {
    input {
        String sample_id
        String partitioning_id
        File inputBAM
        String docker
    }
    
    Int disksize = ceil(4 * size(inputBAM, "GB"))
    
    command <<<
        set -ex
        
        # Create output directory
        mkdir -p partitioned_bams
        cd partitioned_bams/
        
        # Run the partitioning command with error handling
        (
            partition_bam_by_cell_cluster.py \
                --bam ~{inputBAM} \
                --cell_clusters ~{partitioning_id} \
                --output_prefix ~{sample_id} \
                > command_output.log 2>&1
        ) || {
            echo "ERROR: partition_bam_by_cell_cluster.py failed with exit code $?" >&2
            echo "Last 100 lines of output:" >&2
            tail -n 100 command_output.log >&2
            echo "Full command output:" >&2
            cat command_output.log >&2
            exit 1
        }
        
        # List generated BAM files
        echo "Generated BAM files:"
        ls -la *.bam || echo "No BAM files found"
        
        # Verify at least one BAM file was created
        if [ $(ls -1 *.bam 2>/dev/null | wc -l) -eq 0 ]; then
            echo "ERROR: No BAM files were generated" >&2
            exit 1
        fi
    >>>
    
    output {
        Array[File] partitioned_bams = glob("partitioned_bams/*.bam")
    }
    
    runtime {
        docker: docker
        cpu: 1
        memory: "8 GiB"
        disks: "local-disk ~{disksize} HDD"
        maxRetries: 1
    }
    
    meta {
        description: "Partition a BAM file by cell clusters"
    }
    
    parameter_meta {
        sample_id: "Sample identifier used as output prefix"
        partitioning_id: "Partitioning identifier for cell cluster assignments"
        inputBAM: "Input BAM file to partition"
        docker: "Docker image with partition_bam_by_cell_cluster.py"
        partitioned_bams: "Output BAM files, one per cluster"
    }
}
