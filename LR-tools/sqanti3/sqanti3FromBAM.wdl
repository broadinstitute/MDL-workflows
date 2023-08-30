version 1.0


task splitBAMPerChromosomeTask {
    input {
        File inputBAM
        File inputBAMIndex
        Int memoryGB
        String outputType
        String docker
    }

    command <<<
        mkdir -p split_dir
        split_bam_per_chromosome.sh ~{inputBAM} split_dir ~{outputType}
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_dir/*.*am")
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*2 + 10) + " HDD"
        docker: docker
    }
}


task backformatBAMTask {
    input {
        File inputBAM
        String outputType = "sam"       # sam or bam
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    command <<<
        baseBamName=$(basename ~{inputBAM} | sed 's/\(.*\)\..*/\1/')
        # echo "${baseBamName}.~{outputType}" > output_name.txt

        reformat.sh \
            in=~{inputBAM} \
            out=${baseBamName}.backformatted.~{outputType} \
            sam=1.3
    >>>

    output {
        File backformatedBAM = select_first(glob("*.backformatted.~{outputType}"))
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*10 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_CTATLRTask {
    input {
        File inputSAM
        Int memoryGB
        Boolean allowNonPrimary
        # Int diskSizeGB
        String docker
    }

    String alignmentGTF_name = basename("~{inputSAM}", ".sam")
    String extra_arg = if allowNonPrimary then "--allow_non_primary" else ""

    command <<<
        baseSamName=$(basename ~{inputSAM} | sed 's/\(.*\)\..*/\1/')

        SAM_to_gxf.pl --format gtf ~{extra_arg} \
            --sam ~{inputSAM} \
            > temp.sam_to_gxf
        grep -P "[A-z0-1]" temp.sam_to_gxf > {alignmentGTF_name}.gtf
    >>>

    output {
        File alignmentGTF = "~{alignmentGTF_name}"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputSAM, "GB")*3 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_cDNACupcakeTask {
    input {
        File inputSAM
        File referenceFasta
        Boolean correctFasta = false
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    String extra_arg = if correctFasta then "--fasta_correction" else ""
    String alignmentGTF_name = basename("~{inputSAM}", ".sam")

    command <<<
        # baseSamName=$(basename ~{inputSAM} | sed 's/\(.*\)\..*/\1/')

        convert_SAM_to_GTF_for_SQANTI3.py \
            --sam_file ~{inputSAM} \
            --output_prefix ~{alignmentGTF_name} \
            --reference_genome ~{referenceFasta} ~{extra_arg}
    >>>

    output {
        File alignmentGTF = "~{alignmentGTF_name}.gtf"
        File? correctedFasta = "~{alignmentGTF_name}.corrected.fasta"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputSAM, "GB")*3 + 10) + " HDD"
        docker: docker
    }
}


task concatenateGTFsTask {
    input {
        String sampleName
        Array[File] files
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    command <<<
        /scripts/concate_gtfs_and_tag_duplicates.py -o ~{sampleName}.gtf '~{sep="' '" files}'
    >>>

    output {
        File concatenatedGTF = "~{sampleName}.gtf"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk 100 HDD"
        docker: docker
    }
}


task sqantiTask {
    input {
        File inputGTF
        File referenceGTF
        File referenceFasta
        File cagePeak
        File polyAMotifs
        Int cpu
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }


    command <<<
        bash ~{monitoringScript} > monitoring.log &

        sqanti3_qc.py \
            --report both \
            --chunks ~{cpu} \
            --dir sqanti_out_dir \
            --CAGE_peak ~{cagePeak} \
            --polyA_motif_list ~{polyAMotifs} \
            --skipORF \
            --window 20 \
            --isoform_hits \
            ~{inputGTF} \
            ~{referenceGTF} \
            ~{referenceFasta}

            find sqanti_out_dir/ -maxdepth 1 -type f ! -name "*.pdf" -exec gzip {} +
    >>>

    output {
        Array[File] sqantiOutputs = glob("sqanti_out_dir/*")
        File sqantiClassificationTSV = select_first(glob("sqanti_out_dir/*_classification.txt.gz"))
        File sqantiJunctionsTSV = select_first(glob("sqanti_out_dir/*_junctions.txt.gz"))
        File sqantiReportPDF = select_first(glob("sqanti_out_dir/*_SQANTI3_report.pdf"))
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow sqanti3FromBam {

    meta {
        description: "Run Sqanti3 classification of non-assembled reads on an input BAM with the alignments without having to rerun the alignment."
    }

    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        String conversionMethod = "cDNACupcake"
        File referenceGTF
        File referenceFasta
        File cagePeak
        File polyAMotifs
        Boolean allowNonPrimary = true
        Int cpu = 8
        Int memoryGB = 128
        Int diskSizeGB = 500
    }

    String docker = "us-east4-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus@sha256:0da748835f3b95056aa4be3a831c57950a8587c67c17b214a95a53cc94dc3805"

    String outputType = if (conversionMethod == "cDNACupcake") then "sam" else "bam" 
    call splitBAMPerChromosomeTask {
        input:
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            outputType = outputType,
            memoryGB = memoryGB,
            docker = docker
    }

    
    scatter(chromosomeBAM in splitBAMPerChromosomeTask.chromosomeBAMs) {
        if (conversionMethod == "CTAT-LR") {
            call backformatBAMTask {
                input:
                    inputBAM = chromosomeBAM,
                    memoryGB = memoryGB,
                    docker = docker
            }

            call convertSAMtoGTF_CTATLRTask {
                input:
                    inputSAM = backformatBAMTask.backformatedBAM,
                    memoryGB = memoryGB,
                    allowNonPrimary = allowNonPrimary,
                    docker = docker
            }
        }

        if (conversionMethod == "cDNACupcake") {
            call convertSAMtoGTF_cDNACupcakeTask {
                input:
                    inputSAM = chromosomeBAM,
                    referenceFasta = referenceFasta,
                    memoryGB = memoryGB,
                    docker = docker
            }
        }

        File convertedGTF = select_first([convertSAMtoGTF_CTATLRTask.alignmentGTF, convertSAMtoGTF_cDNACupcakeTask.alignmentGTF])
    }

    call concatenateGTFsTask {
        input:
            sampleName = sampleName,
            files = convertedGTF,
            memoryGB = 8,
            docker = docker
    }

    call sqantiTask {
        input:
            inputGTF = concatenateGTFsTask.concatenatedGTF,
            referenceGTF = referenceGTF,
            referenceFasta = referenceFasta,
            cagePeak = cagePeak,
            polyAMotifs = polyAMotifs,
            cpu = cpu,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB,
            docker = docker
    }

    output {
        Array[File] sqantiOutputs = sqantiTask.sqantiOutputs
        File sqantiClassificationTSV = sqantiTask.sqantiClassificationTSV
        File sqantiJunctionsTSV = sqantiTask.sqantiJunctionsTSV
        File sqantiReportPDF = sqantiTask.sqantiReportPDF
    }
}

