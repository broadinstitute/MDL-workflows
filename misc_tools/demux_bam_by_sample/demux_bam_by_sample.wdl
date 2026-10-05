version 1.0

workflow Demux_Bam_By_Sample {
    input {
        File   bam
        File   metadata_csv
        String pool_name
        String barcode_tag = "CB"
        Boolean trim_barcode_suffix = true
        String barcode_column = "barcode"
        String donor_column   = "donor_id"
        String sample_column  = "sample"
        String pool_column    = "orig.ident"
        Int    threads = 4
        Int    preemptible = 3
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/demux_bam_by_sample:0.1.0"
    }

    call Demux {
        input:
            bam = bam, metadata_csv = metadata_csv, pool_name = pool_name,
            barcode_tag = barcode_tag, trim_barcode_suffix = trim_barcode_suffix,
            barcode_column = barcode_column, donor_column = donor_column,
            sample_column = sample_column, pool_column = pool_column,
            threads = threads,
            preemptible = preemptible, docker = docker
    }

    output {
        # All four arrays share the same order (sorted by BAM file name).
        Array[File]   sample_bams    = Demux.sample_bams
        Array[File]   sample_bais    = Demux.sample_bais
        Array[String] donor_ids      = Demux.donor_ids
        Array[String] samples        = Demux.samples
        File          demux_summary  = Demux.demux_summary
        File          demux_totals   = Demux.demux_totals
        File          demux_log      = Demux.demux_log
    }
}

task Demux {
    input {
        File   bam
        File   metadata_csv
        String pool_name
        String barcode_tag
        Boolean trim_barcode_suffix
        String barcode_column
        String donor_column
        String sample_column
        String pool_column
        Int    threads
        Int    preemptible
        String docker
    }

    String prefix = basename(bam, ".bam")
    String trim_arg = if trim_barcode_suffix then "--trim-suffix" else "--no-trim-suffix"
    # Outputs hold the same reads as the input, split in n files.
    Int disk_gb = ceil(size(bam, "GB") * 2.3 + size(metadata_csv, "GB") + 20)

    command <<<
        set -euo pipefail
        mkdir out
        demux_bam_by_sample \
            --bam ~{bam} \
            --metadata ~{metadata_csv} \
            --pool ~{pool_name} \
            --barcode-tag ~{barcode_tag} \
            ~{trim_arg} \
            --barcode-column ~{barcode_column} \
            --donor-column ~{donor_column} \
            --sample-column ~{sample_column} \
            --pool-column ~{pool_column} \
            --threads ~{threads} \
            --output-prefix ~{prefix} \
            --output-dir out 2> demux.log || { cat demux.log >&2; exit 1; }
        cat demux.log >&2

        # Columns of the summary: donor_id, sample, bam, bai.
        # sample_bams / sample_bais are globbed (read_lines would only give path strings, which Cromwell does not delocalize).
        # glob() returns files sorted by name, so donor_ids / samples are written in the same (C-locale, by BAM name) order.
        tail -n +2 out/~{prefix}.demux_summary.tsv | LC_ALL=C sort -t "$(printf '\t')" -k3,3 > summary_by_bam.tsv
        cut -f1 summary_by_bam.tsv > donor_ids.txt
        cut -f2 summary_by_bam.tsv > samples.txt
    >>>

    output {
        Array[File]   sample_bams   = glob("out/*.bam")
        Array[File]   sample_bais   = glob("out/*.bam.bai")
        Array[String] donor_ids     = read_lines("donor_ids.txt")
        Array[String] samples       = read_lines("samples.txt")
        File          demux_summary = "out/~{prefix}.demux_summary.tsv"
        File          demux_totals  = "out/~{prefix}.demux_totals.tsv"
        File          demux_log     = "demux.log"
    }

    runtime {
        docker: docker
        predefinedMachineType: "c3d-highcpu-4"
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }
}
