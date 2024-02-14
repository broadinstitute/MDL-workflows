version 1.0

task LongRNAqcPlottingTask {
    input {
        Array[String] sampleName
        Array[File] classificationFile
        Array[File] junctionFile
        String outputPrefix
        Boolean includeSaturation
        Int maxRetries
    }

    # Calculate total memory required
    Int total_file_size = ceil(size(classificationFile, "GiB") + size(junctionFile, "GiB") + 8)
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
        disks: "local-disk " + total_file_size*2 + " HDD"
        cpu: 1
        memory: total_file_size*2 + " GiB"
        preemptible: maxRetries
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
        Int maxRetries = 3
    }

    call LongRNAqcPlottingTask {
        input:
            sampleName = sampleName,
            classificationFile = classificationFile,
            junctionFile = junctionFile,
            outputPrefix = outputPrefix,
            includeSaturation = includeSaturation,
            maxRetries = maxRetries
    }

    output {
        File QC_plots = LongRNAqcPlottingTask.plots
    }
}

