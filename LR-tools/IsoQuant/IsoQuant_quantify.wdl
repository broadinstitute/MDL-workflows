version 1.0


task run_isoquant_quantify {
    input {
        File inputBAM
        File referenceFasta
        File geneDB
        String dataType
        Boolean noModelConstruction
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-isoquant/lrtools-isoquant-plus@sha256:afad1eba2743f09cc8bddf6f38b99f3b8fd104c67dddccd830ecb2e43ec3deab"
        File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"
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
            --threads ~{numThreads} ~{extra_args}\
            -o isoquant_out
    >>>

    output {
        Array[File] isoquantOutputs = glob("isoquant_out/*")
        File readAssignmentsTsv = select_first(glob("isoquant_out/*.read_assignments.tsv"))
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow isoquant_quantify {
    meta {
        description: "Run IsoQuant quantification (on an already gffutils preprocessed reference geneDB ideally)."
    }

    input {
        File inputBAM
        File referenceFasta
        File geneDB
        String dataType = "pacbio_ccs"
        Boolean noModelConstruction = false
    }

    call run_isoquant_quantify {
        input:
            inputBAM = inputBAM,
            geneDB = geneDB,
            referenceFasta = referenceFasta,
            dataType = dataType
    }

    output {
        Array[File] isoquantOutputs = run_isoquant_quantify.isoquantOutputs
        File readAssignmentsTsv = run_isoquant_quantify.readAssignmentsTsv
        File monitoringLog = run_isoquant_quantify.monitoringLog
    }
}