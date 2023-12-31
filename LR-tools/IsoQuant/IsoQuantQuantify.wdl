version 1.0


task isoquantQuantifyTask {
    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        File referenceFasta
        File geneDB
        String dataType
        Boolean noModelConstruction
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-isoquant/lrtools-isoquant-plus@sha256:bbad9d6cb47bcaa6de76c04d425bd3815d7f4b12f5679dac2eb894aa4ee3f81f"
        Int preemptible_tries
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    # String file_name = basename("~{inputBAM}", ".bam")
    String extra_args = if noModelConstruction then "--no_model_construction" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.3.1/isoquant.py \
            --reference ~{referenceFasta} \
            --genedb ~{geneDB} \
            --bam ~{inputBAM} \
            --data_type ~{dataType} \
            --stranded forward \
            --transcript_quantification all \
            --gene_quantification all \
            --threads ~{numThreads} ~{extra_args} \
            --labels ~{sampleName} \
            --prefix ~{sampleName} \
            -o isoquant_output

            find isoquant_output/~{sampleName}/ -maxdepth 1 -type f -exec gzip {} +
    >>>

    output {
        Array[File] isoquantOutputs = glob("isoquant_output/~{sampleName}/*.gz")
        File readAssignmentsTSV = select_first(glob("isoquant_output/~{sampleName}/*.read_assignments.tsv.gz"))
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
        preemptible: preemptible_tries
    }
}


workflow isoquantQuantify {
    meta {
        description: "Run IsoQuant quantification (on an already gffutils preprocessed reference geneDB ideally)."
    }

    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        File referenceFasta
        File geneDB
        String dataType = "pacbio_ccs"
        Boolean noModelConstruction = false
        Int preemptible_tries = 3
    }

    call isoquantQuantifyTask {
        input:
            sampleName = sampleName,
            inputBAM = inputBAM,
            inputBAMIndex = inputBAMIndex,
            geneDB = geneDB,
            referenceFasta = referenceFasta,
            dataType = dataType,
            noModelConstruction = noModelConstruction,
            preemptible_tries = preemptible_tries
    }

    output {
        Array[File] isoquantOutputs = isoquantQuantifyTask.isoquantOutputs
        File readAssignmentsTSV = isoquantQuantifyTask.readAssignmentsTSV
        # File monitoringLog = isoquantQuantifyTask.monitoringLog
    }
}