version 1.0


task run_liqa_quantify {
    input {
        File inputBAM
        File prepared_reference
        Int max_distance = 20 # The maximum length of an alignment error at exon boundary. Recommend: 20.
        Int f_weight = 1 # The weight for bias correction in isoform usage estimation. Recommend: 1
        Int cpu = 16
        Int numThreads = 32
        Int memoryGB = 128
        Int diskSizeGB = 500
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lrtools-liqa/lrtools-liqa@sha256:c0fb62e176893bb30d63c43bd9d108a3d6e3c8a24c3ba5c0110b4d4ff2b73476"
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }

    String file_name = basename("~{inputBAM}", ".bam")

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        liqa -task quantify \
        -refgene ~{prepared_reference} \
        -bam ~{inputBAM} \
        -max_distance ~{max_distance} \
        -f_weight ~{f_weight} \
        -out ~{file_name}.liqa_out.tsv
    >>>

    output {
        File liqaOutput = "~{file_name}.liqa_out.tsv"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow liqa_quantify {
    meta {
        description: "Run LIQA quantification on an already preprocessed reference."
    }

    input {
        File inputBAM
        File prepared_reference
        Int max_distance = 20 # The maximum length of an alignment error at exon boundary. Recommend: 20.
        Int f_weight = 1 # The weight for bias correction in isoform usage estimation. Recommend: 1
    }

    call run_liqa_quantify {
        input:
            inputBAM = inputBAM,
            prepared_reference = prepared_reference,
            max_distance = max_distance,
            f_weight = f_weight
    }

    output {
        File liqaOutput = run_liqa_quantify.liqaOutput
        File monitoringLog = run_liqa_quantify.monitoringLog
    }

}