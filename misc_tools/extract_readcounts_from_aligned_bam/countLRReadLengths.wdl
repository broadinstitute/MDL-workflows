version 1.0

## WDL wrapper for counting read lengths from BAM/CRAM files
## This workflow processes mapped reads only using parallel per-contig processing
## Author: cgeorges@broadinstitute.org (original script)

workflow CountReadLengths {
    input {
        File input_bam
        File? input_bam_index
        File python_script
        String output_basename = "read_length_counts"
        Int workers = 4
    }

    call CountReadLengthsTask {
        input:
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            python_script = python_script,
            output_basename = output_basename,
            workers = workers
    }

    output {
        File read_length_counts = CountReadLengthsTask.read_length_counts
    }
}

task CountReadLengthsTask {
    input {
        File input_bam
        File? input_bam_index
        File python_script
        String output_basename
        Int workers
        
        # Runtime parameters
        Int cpu = 4
        Int memory_gb = 8
        Int disk_gb = 100
        String docker = "quay.io/biocontainers/pysam:0.21.0--py39h2bbff1b_1"
    }

    String output_filename = "${output_basename}.pkl"
    
    command <<<
        set -euo pipefail
        
        # Copy the input Python script to the working directory
        cp "~{python_script}" count_read_lengths.py
        
        # Make the script executable
        chmod +x count_read_lengths.py
        
        # Run the script with mapped reads only (no --include-unmapped flag)
        python3 count_read_lengths.py \
            --bam "~{input_bam}" \
            --out "~{output_filename}" \
            --workers ~{workers}
    >>>

    output {
        File read_length_counts = output_filename
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${memory_gb} GB"
        disks: "local-disk ${disk_gb} HDD"
    }
}
