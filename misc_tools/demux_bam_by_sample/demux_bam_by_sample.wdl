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
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/demux_bam_by_sample:0.2.0"
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
        # All arrays share the same order (sorted by BAM file name).
        Array[File]   demux_sample_bams    = Demux.sample_bams
        Array[File]   demux_sample_bais    = Demux.sample_bais
        Array[String] demux_donor_ids      = Demux.donor_ids
        Array[String] demux_samples        = Demux.samples
        Array[Int]    demux_reads    = Demux.demux_reads
        Array[Int]    demux_cells    = Demux.demux_cells
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
    # Outputs hold the same reads as the input, split in n files; the factor 3 leaves room to sort one output when the input is unsorted.
    Int disk_gb = ceil(size(bam, "GB") * 3 + size(metadata_csv, "GB") + 20)

    command <<<
        set -euo pipefail
        # Declared as a bash function so the same call can be repeated with --no-index.
        demux() {
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
                --output-dir out "$@"
        }

        # Try with indexes first (needs coordinate order; a header that does not say SO:coordinate is accepted and checked read by read).
        # Exit code 3 = a read broke the order: redo without indexes, then sort and index each sample BAM here.
        mkdir out
        set +e
        demux 2> demux.log
        rc=$?
        set -e
        if [ "$rc" -eq 3 ]; then
            echo "[info] input is not coordinate-sorted: rerunning with --no-index, then sorting and indexing the outputs" >> demux.log
            mv demux.log demux.attempt1.log
            rm -rf out && mkdir out
            demux --no-index 2> demux.log || { cat demux.attempt1.log demux.log >&2; exit 1; }
            for b in out/*.bam; do
                samtools sort --no-PG --write-index -@ ~{threads} -T "$b.tmp" -o "$b.sorted##idx##$b.sorted.bai" "$b"
                mv "$b.sorted" "$b"
                mv "$b.sorted.bai" "$b.bai"
            done
            cat demux.attempt1.log demux.log > demux.both.log && mv demux.both.log demux.log && rm -f demux.attempt1.log
        elif [ "$rc" -ne 0 ]; then
            cat demux.log >&2
            exit "$rc"
        fi
        cat demux.log >&2

        # Columns of the summary: donor_id, sample, bam, bai, n_barcodes, n_reads, status, n_barcodes_seen.
        # sample_bams / sample_bais are globbed (read_lines would only give path strings, which Cromwell does not delocalize).
        # glob() returns files sorted by name, so donor_ids / samples are written in the same (C-locale, by BAM name) order.
        tail -n +2 out/~{prefix}.demux_summary.tsv | LC_ALL=C sort -t "$(printf '\t')" -k3,3 > summary_by_bam.tsv
        cut -f1 summary_by_bam.tsv > donor_ids.txt
        cut -f2 summary_by_bam.tsv > samples.txt
        # Reads assigned to each sample, and the number of its barcodes (cells) actually observed in the reads.
        cut -f6 summary_by_bam.tsv > demux_reads.txt
        cut -f8 summary_by_bam.tsv > demux_cells.txt
    >>>

    output {
        Array[File]   sample_bams   = glob("out/*.bam")
        Array[File]   sample_bais   = glob("out/*.bam.bai")
        Array[String] donor_ids     = read_lines("donor_ids.txt")
        Array[String] samples       = read_lines("samples.txt")
        Array[Int]    demux_reads   = read_lines("demux_reads.txt")
        Array[Int]    demux_cells   = read_lines("demux_cells.txt")
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
