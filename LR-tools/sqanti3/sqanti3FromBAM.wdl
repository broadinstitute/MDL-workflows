version 1.0


task splitGTFPerChromosomeTask {
    input {
        File inputGTF
        String chromosomesList
        Int memoryGB = 16
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &
        
        mkdir -p split_dir
        split_gtf_per_chromosome.sh ~{inputGTF} split_dir ~{chromosomesList}
    >>>

    output {
        Array[File] chromosomeGTFs = glob("split_dir/*.gtf")
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputGTF, "GB")*2 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_CTATLRTask {
    input {
        File inputBAM
        File inputBAMIndex
        Int memoryGB = 16
        Boolean allowNonPrimary
        # Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    String baseBamName = basename("~{inputBAM}", ".bam")
    String extra_arg = if allowNonPrimary then "--allow_non_primary" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        reformat.sh \
            in=~{inputBAM} \
            out=${baseBamName}.backformatted.sam \
            sam=1.3

        SAM_to_gxf.pl --format gtf ~{extra_arg} \
            --sam ${baseBamName}.backformatted.sam \
            > temp.sam_to_gxf
        grep -P "[A-z0-1]" temp.sam_to_gxf > {baseBamName}.gtf
    >>>

    output {
        File alignmentGTF = "~{baseBamName}"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*20 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_cDNACupcakeTask {
    input {
        File inputBAM
        File inputBAMIndex
        File referenceFasta
        Boolean correctFasta = false
        Boolean allowNonPrimary = true
        Int memoryGB = 32
        # Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    String extra_arg = if allowNonPrimary then "--allow_non_primary" else ""
    String extra_arg2 = if correctFasta then "--fasta_correction" else ""
    String baseBamName = basename("~{inputBAM}", ".bam")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        samtools view ~{inputBAM} > tmp.sam

        convert_SAM_to_GTF_for_SQANTI3.py \
            --sam_file  tmp.sam \
            --output_prefix ~{baseBamName} \
            --reference_genome ~{referenceFasta} ~{extra_arg} ~{extra_arg2}

    >>>

    output {
        File alignmentGTF = "~{baseBamName}.gtf"
        File? correctedFasta = "~{baseBamName}.corrected.fasta"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*30 + 10) + " SSD"
        docker: docker
    }
}


task concatenateSqantiOutputsTask {
    input {
        String sampleName
        Array[File] classificationFiles
        Array[File] junctionFiles
        Array[File] correctedFastaFiles
        Array[File] correctedGTFFiles
        Int memoryGB = 16
        # Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        cat '~{sep="' '" classificationFiles}' | gzip > ~{sampleName}_classification.tsv.gz
        cat '~{sep="' '" junctionFiles}' | gzip > ~{sampleName}_junctions.tsv.gz
        cat '~{sep="' '" correctedFastaFiles}' | gzip > ~{sampleName}_corrected.fasta.gz
        cat '~{sep="' '" correctedGTFFiles}' | gzip > ~{sampleName}_corrected.gtf.gz
    >>>

    output {
        File concatenatedClassification = "~{sampleName}_classification.tsv.gz"
        File concatenatedJunctions = "~{sampleName}_junctions.tsv.gz"
        File concatenatedCorrectedFasta = "~{sampleName}_corrected.fasta.gz"
        File concatenatedCorrectedGTF = "~{sampleName}_corrected.gtf.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk 500 HDD"
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
        Int memoryGB
        Int diskSizeGB
        String docker
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }


    command <<<
        bash ~{monitoringScript} > monitoring.log &

        sqanti3_qc.py \
            --report skip \
            --dir sqanti_out_dir \
            --CAGE_peak ~{cagePeak} \
            --polyA_motif_list ~{polyAMotifs} \
            --skipORF \
            --window 20 \
            --isoform_hits \
            ~{inputGTF} \
            ~{referenceGTF} \
            ~{referenceFasta}

            # find sqanti_out_dir/ -maxdepth 1 -type f ! -name "*.pdf" -exec gzip {} +
    >>>

    output {
        Array[File] sqantiOutputs = glob("sqanti_out_dir/*")
        File sqantiClassificationTSV = select_first(glob("sqanti_out_dir/*_classification.txt"))
        File sqantiJunctionsTSV = select_first(glob("sqanti_out_dir/*_junctions.txt"))
        File sqantiCorrectedFasta = select_first(glob("sqanti_out_dir/*_corrected.fasta"))
        File sqantiCorrectedGTF = select_first(glob("sqanti_out_dir/*_corrected.gtf"))
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: 1
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
        String chromosomesList # comma seprarated
        String conversionMethod = "cDNACupcake"
        File referenceGTF
        File referenceFasta
        File cagePeak
        File polyAMotifs
        Boolean allowNonPrimary = true
        Int memoryGB = 32
        Int diskSizeGB = 256
    }

    String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus@sha256:a7f116572bc67f5c80165ab90b9baea0cb284d3a0b1b7196601784c781108d92"

    if (conversionMethod == "CTAT-LR") {
        call convertSAMtoGTF_CTATLRTask {
            input:
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                # memoryGB = memoryGB,
                allowNonPrimary = allowNonPrimary,
                docker = docker
        }
    }

    if (conversionMethod == "cDNACupcake") {
        call convertSAMtoGTF_cDNACupcakeTask {
            input:
                inputBAM = inputBAM,
                inputBAMIndex = inputBAMIndex,
                referenceFasta = referenceFasta,
                # memoryGB = memoryGB,
                docker = docker
        }
    }

    File convertedGTF = select_first([convertSAMtoGTF_CTATLRTask.alignmentGTF, convertSAMtoGTF_cDNACupcakeTask.alignmentGTF])
    
    call splitGTFPerChromosomeTask {
        input:
            inputGTF = convertedGTF,
            chromosomesList = chromosomesList,
            # memoryGB = memoryGB,
            docker = docker
    }


    scatter(chromosomeGTF in splitGTFPerChromosomeTask.chromosomeGTFs) {
        call sqantiTask {
            input:
                inputGTF = chromosomeGTF,
                referenceGTF = referenceGTF,
                referenceFasta = referenceFasta,
                cagePeak = cagePeak,
                polyAMotifs = polyAMotifs,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker
        }
    }


    call concatenateSqantiOutputsTask {
        input:
            sampleName = sampleName,
            classificationFiles = sqantiTask.sqantiClassificationTSV,
            junctionFiles = sqantiTask.sqantiJunctionsTSV,
            correctedFastaFiles = sqantiTask.sqantiCorrectedFasta,
            correctedGTFFiles = sqantiTask.sqantiCorrectedGTF,
            memoryGB = 8,
            docker = docker
    }


    output {
        File sqantiClassificationTSV = concatenateSqantiOutputsTask.concatenatedClassification
        File sqantiJunctionsTSV = concatenateSqantiOutputsTask.concatenatedJunctions
        File sqantiCorrectedFasta = concatenateSqantiOutputsTask.concatenatedCorrectedFasta
        File sqantiCorrectedGTF = concatenateSqantiOutputsTask.concatenatedCorrectedGTF
    }
}

