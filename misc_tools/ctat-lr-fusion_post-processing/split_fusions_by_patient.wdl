version 1.0

workflow SplitFusionsByPatient {
    input {
        File input_tsv
        File path_split_fusions_by_id_py
        String base_output_dir = "fusions_by_pool"
    }

    call SplitFusions {
        input:
            input_tsv = input_tsv,
            path_split_fusions_by_id_py = path_split_fusions_by_id_py,
            base_output_dir = base_output_dir
    }

    output {
        Directory split_files = SplitFusions.output_directory
        Array[File] split_tsvs = SplitFusions.split_tsv_files
    }
}

task SplitFusions {
    input {
        File input_tsv
        File path_split_fusions_by_id_py
        String base_output_dir
    }

    command <<<
        # Install pandas
        pip install pandas
        
        # Create the base output directory
        mkdir -p ~{base_output_dir}

        # Copy the script to working directory
        cp ~{path_split_fusions_by_id_py} split_fusions_by_patient.py
        
        # Run the Python script
        python3 ./split_fusions_by_patient.py ~{input_tsv} ~{base_output_dir}
        
        # Find all generated TSV files for output
        find ~{base_output_dir} -name "*.tsv" > tsv_files.txt
    >>>

    output {
        Directory output_directory = base_output_dir
        Array[File] split_tsv_files = glob("${base_output_dir}/*/*.tsv")
    }

    runtime {
        docker: "python:3.9-slim"
        memory: "8 GB"
        cpu: 2
        disks: "local-disk 32 SSD"
    }

    meta {
        description: "Split fusion calls TSV by patient assignment"
    }
}
