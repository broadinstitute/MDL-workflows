version 1.0

task PlotMatchRatesAndIndelRates {
    input {
        Array[File] bamFiles
        Array[File] bamIndexes
        Array[String] sampleNames
        Int ?maxRetries = 3
    }

   Int total_file_size = ceil(size(bamFiles, "GiB"))

    command <<<
        violin_plot_indel_rates.py "~{sep=',' sampleNames}" "~{sep=',' bamFiles}"
    >>>

    runtime {
        docker: "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-qc/longrnaqc_readaccuracyanalyzer"
        disks: "local-disk " + ceil(total_file_size*1.1 + 5) + " SSD"
        memory: (total_file_size/2 + 4) + " GiB"
        cpu: 1
        preemptible: maxRetries
    }

    output {
        File violinPlots = "violin_plots.png"
        File phredViolinPlots = "phred_violin_plots.png"
    }
}

workflow LongRNAqc_readAccuracyAnalyzer {
    meta {
        description: "Generate violin plots of in/del rates, phred scores, and proportion of sequence from BAMs."
    }

    input {
        Array[File] bamFiles
        Array[File] bamIndexes
        Array[String] sampleNames
        Int ?maxRetries = 3
    }

    call PlotMatchRatesAndIndelRates {
        input:
            bamFiles = bamFiles,
            bamIndexes = bamIndexes,
            sampleNames = sampleNames,
            maxRetries = maxRetries
    }

    output {
        File violinPlots = PlotMatchRatesAndIndelRates.violinPlots
        File phredViolinPlots = PlotMatchRatesAndIndelRates.phredViolinPlots
    }
}