version 1.0

# Pre-alignment collapse (step 2): reads sharing CB+UMI are collapsed into consensus records before alignment.
# Input = the bucket-ordered unaligned BAMs of step 1 (BarcodeExtractAndCorrect with out_format = "bam" and buckets > 0) and their
# .bkt indexes. A planning task cuts the bucket range into contiguous intervals of about target_reads_per_shard reads (summed over
# every input file, so several runs of one library can be collapsed together if their step 1 used the same buckets); each shard
# reads its byte range of every BAM straight from GCS (no localization) and writes cDNA-oriented consensus records (FASTQ, BGZF)
# for the alignment step. Optionally the first shard is first run as a calibration that measures the error profile on its own
# families and writes the tie tables every shard then uses.

workflow PreAlignmentCollapse {
    input {
        # step-1 outputs: fastq_pairs (bam) and bucket_indexes, same nesting and order
        Array[Array[File]] step1_bams
        Array[Array[File]] step1_bucket_indexes
        String sample_name
        String priming                       # 5p | 3p (library); reads are processed as sequenced, records are written in cDNA orientation

        Int target_reads_per_shard = 90000000 # ~0.55 GB of collapse memory per million reads: 90 M fits a 64 GB machine
        Boolean calibrate = true             # tie tables from the first shard's families (else the tables compiled into collapse)
        Int ref_end_tol = 2
        String extra_collapse_args = ""

        String machine_type = "c3d-standard-16"
        Int threads = 16
        Int preemptible = 3
        String docker
    }

    Array[File] bams = flatten(step1_bams)
    Array[File] bucket_indexes = flatten(step1_bucket_indexes)

    call Plan_Shards {
        input:
            bams = bams,
            bucket_indexes = bucket_indexes,
            target_reads = target_reads_per_shard,
            docker = docker
    }

    if (calibrate) {
        call Collapse_Shard as calibration {
            input:
                manifest = Plan_Shards.manifests[0],
                bams = bams,
                priming = priming,
                out_name = sample_name + ".calibration",
                write_tie_tables = true,
                ref_end_tol = ref_end_tol,
                extra_args = extra_collapse_args,
                machine_type = machine_type,
                threads = threads,
                preemptible = preemptible,
                docker = docker
        }
    }

    scatter (i in range(length(Plan_Shards.manifests))) {
        call Collapse_Shard {
            input:
                manifest = Plan_Shards.manifests[i],
                bams = bams,
                priming = priming,
                out_name = sample_name + ".shard" + (i + 1),
                tie_tables_in = calibration.tie_tables,
                ref_end_tol = ref_end_tol,
                extra_args = extra_collapse_args,
                machine_type = machine_type,
                threads = threads,
                preemptible = preemptible,
                docker = docker
        }
    }

    output {
        Array[File] collapsed_fastqs = Collapse_Shard.records
        Array[File] collapse_stats = Collapse_Shard.log
        File shard_plan = Plan_Shards.summary
        File? tie_tables = calibration.tie_tables
        File? error_profile = calibration.profile
    }
}

task Plan_Shards {
    input {
        Array[File] bams           # not localized: only their paths go into the manifests
        Array[File] bucket_indexes # small, localized
        Int target_reads
        String docker
    }

    parameter_meta {
        bams: { localization_optional: true }
    }

    command <<<
        set -euo pipefail
        paste ~{write_lines(bams)} ~{write_lines(bucket_indexes)} > files.tsv
        plan_buckets.py --files files.tsv --target-reads ~{target_reads} --prefix shard_
        ls shard_[0-9][0-9][0-9][0-9].tsv > manifests.list
    >>>

    output {
        Array[File] manifests = read_lines("manifests.list")
        File summary = "shard_summary.tsv"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
    }
}

task Collapse_Shard {
    input {
        File manifest
        Array[File] bams # read by byte range from GCS (paths in the manifest); declared so the workflow tracks them as files
        String priming
        String out_name
        File? tie_tables_in
        Boolean write_tie_tables = false
        Int ref_end_tol
        String extra_args
        String machine_type
        Int threads
        Int preemptible
        String docker
    }

    parameter_meta {
        bams: { localization_optional: true }
    }

    command <<<
        set -euo pipefail
        if grep -q $'^gs://' ~{manifest}; then
          # short-lived token of the VM's service account (metadata server; no key); htslib sends it with each gs:// request
          export GCS_OAUTH_TOKEN=$(curl -sf -H 'Metadata-Flavor: Google' \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json, sys; print(json.load(sys.stdin)["access_token"])')
        fi
        if [[ "~{write_tie_tables}" == true ]]; then
          collapse -p ~{priming} -t ~{threads} --ref-end-tol ~{ref_end_tol} ~{extra_args} --ranges ~{manifest} \
            --write-tie-tables ~{out_name}.tie_tables.txt --profile-out ~{out_name}.profile.tsv -o /dev/null 2> ~{out_name}.log
          : > ~{out_name}.fastq.gz
        else
          collapse -p ~{priming} -t ~{threads} --ref-end-tol ~{ref_end_tol} ~{extra_args} --ranges ~{manifest} \
            ~{if defined(tie_tables_in) then "--tie-tables " + tie_tables_in else ""} -o ~{out_name}.fastq.gz 2> ~{out_name}.log
        fi
        tail -n 3 ~{out_name}.log
    >>>

    output {
        File records = "~{out_name}.fastq.gz"
        File log = "~{out_name}.log"
        File? tie_tables = "~{out_name}.tie_tables.txt"
        File? profile = "~{out_name}.profile.tsv"
    }

    # cpu/memory omitted: GCP Batch mishandles predefinedMachineType together with explicit cpu/memory
    runtime {
        docker: docker
        predefinedMachineType: machine_type
        disks: "local-disk 50 SSD"
        preemptible: preemptible
    }
}
