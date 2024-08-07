version 1.0

task LongRNAqcPlottingTask {
    input {
        Array[String] sampleName
        Array[File] classificationFile
        Array[File] junctionFile
        String outputPrefix
        Boolean includeSaturation
        Int? memoryGB
        Int preemptible
    }

    # Calculate total memory required
    Int total_file_size = ceil(size(classificationFile, "GiB") + size(junctionFile, "GiB") + 8)
    Int total_classification_file_size = ceil(size(classificationFile, "GiB") + size(junctionFile, "GiB"))

    Int memory_use = select_first([memoryGB, total_classification_file_size*7 + 8])

    String saturation_arg = if includeSaturation then "True" else "False"
    
    command {
        report_multisample_shortform.R \
            '~{sep="," sampleName}' \
            '~{sep="," classificationFile}' \
            '~{sep="," junctionFile}' \
            ~{outputPrefix} \
            ~{saturation_arg}
    }

    output {
        File plots = select_first(glob("~{outputPrefix}*"))
    }

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plotting"
        bootDiskSizeGb: 30
        disks: "local-disk " + total_file_size*2 + " HDD"
        cpu: 1
        memory: memory_use + " GiB"
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
        Array[File] junctionFile
        String outputPrefix
        Boolean includeSaturation
        Int preemptible = 1
    }

    call LongRNAqcPlottingTask {
        input:
            sampleName = sampleName,
            classificationFile = classificationFile,
            junctionFile = junctionFile,
            outputPrefix = outputPrefix,
            includeSaturation = includeSaturation,
            preemptible = preemptible
    }

    output {
        File QC_plots = LongRNAqcPlottingTask.plots
    }
}

