version 1.0


task isoquantMakeGeneDBTask {
    input {
        File gtfToDB
        Boolean isCompleteGeneDB
        Int cpu = 1
        Int memoryGB = 16
        Int diskSizeGB = 50
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-isoquant/lrtools-isoquant-plus@sha256:b1f6d4e6db845372e31d81dfbf8cd5e06f6f68dd1b6cd605a7a2b5f6211f3aa1"
        Int preemptible_tries
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    String extra_args = if isCompleteGeneDB then "--complete_genedb" else ""

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        /usr/local/src/IsoQuant-3.4.1/isoquant_prepare_genedb.py \
            --genedb ~{gtfToDB} \
            --genedb_output ./ \
            ~{extra_args}
    >>>

    output {
        File geneDB = select_first(glob("*.db"))
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


workflow isoquantMakeGeneDB {
    meta {
        description: "Run IsoQuant only to make a .db file from a .gtf file."
    }

    input {
        File gtfToDB
        Boolean isCompleteGeneDB = false
        Int preemptible_tries = 3
    }

    call isoquantMakeGeneDBTask {
        input:
            gtfToDB = gtfToDB,
            isCompleteGeneDB = isCompleteGeneDB,
            preemptible_tries = preemptible_tries
    }

    output {
        File geneDB = isoquantMakeGeneDBTask.geneDB
        # File monitoringLog = isoquantMakeGeneDBTask.monitoringLog
    }
}