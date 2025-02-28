version 1.0

task LongRNAqcPlottingTask {
    input {
        Array[String] sampleName
        Array[File] classificationFile
        String outputPrefix
        Int preemptible
    }

    # Calculate total memory required
    Int total_file_size = ceil(size(classificationFile, "GiB") + 8)

    command {
        /scripts/LongRNAqc_classification_plots.py \
            --sample_names '~{sep="," sampleName}' \
            --classification_files '~{sep="," classificationFile}' \
            --output ~{outputPrefix} \
            --type sqanti
    }

    output {
        File QC_categories_plots = "~{outputPrefix}_categories.pdf"
        File QC_read_lengths_plots = "~{outputPrefix}_read_length_distributions.pdf"
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-qc/lrtools-qc-plotting"
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size*2 + " HDD"
        cpu: 1
        memory: "4 GiB"
        preemptible: preemptible
    }
}


workflow LongRNAqcPlotting {
    meta {
        description: "Generate multi-sample QC plots using Sqanti3 outputs."
    }

    input {
        Array[String] sampleName
        Array[File] classificationFile
        String outputPrefix
        Int preemptible = 1
    }

    call LongRNAqcPlottingTask {
        input:
            sampleName = sampleName,
            classificationFile = classificationFile,
            outputPrefix = outputPrefix,
            preemptible = preemptible
    }

    output {
        File QC_categories_plots = LongRNAqcPlottingTask.QC_categories_plots
        File QC_read_lengths_plots = LongRNAqcPlottingTask.QC_read_lengths_plots
    }
}

