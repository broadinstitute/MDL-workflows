version 1.0

# Cross-sample comparison of the tables written by Bam_Dedup_QC: overlay plots + one summary row per sample and stage.
# Takes plain parallel arrays (one entry per sample, same order in every array), so it does not depend on how the
# arrays are gathered (a set of sets, a sample set, a JSON file ...). Pre and post are independent: give either or both.
workflow Bam_Dedup_QC_Compare {
    input {
        String comparison_name

        Array[String] pre_sample_names = []
        Array[File] pre_aligned_length_tsvs = []
        Array[File] pre_error_by_position_tsvs = []
        Array[File] pre_read_error_rates_tsvs = []  # optional (per-read plot)
        Array[File] pre_cb_counts_tsv_gzs = []
        Array[File] pre_summary_jsons = []

        Array[String] post_sample_names = []
        Array[File] post_aligned_length_tsvs = []
        Array[File] post_dc_tsvs = []
        Array[File] post_error_by_position_tsvs = []
        Array[File] post_read_error_rates_tsvs = []
        Array[File] post_cb_counts_tsv_gzs = []
        Array[File] post_summary_jsons = []

        String error_plot_style = "box"
        Int cell_min_reads = 0  # 0 = estimate each dataset's cell cutoff from its knee
        String cell_cutoff_method = "steepest"
        Float error_axis_max = 0.04  # common upper limit of every error-rate axis (fraction of reference bases)
        Int memory_gb = 16
        Int preemptible = 3
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/bam_dedup_qc:0.1.0"
    }

    # Fail before any file is localized: every table array must have one entry per sample name.
    call CheckInputs {
        input:
            n_pre_names = length(pre_sample_names), n_pre_tables = [length(pre_aligned_length_tsvs), length(pre_error_by_position_tsvs), length(pre_cb_counts_tsv_gzs), length(pre_summary_jsons)],
            n_pre_optional = [length(pre_read_error_rates_tsvs)],
            n_post_names = length(post_sample_names), n_post_tables = [length(post_aligned_length_tsvs), length(post_dc_tsvs), length(post_error_by_position_tsvs), length(post_cb_counts_tsv_gzs), length(post_summary_jsons)],
            n_post_optional = [length(post_read_error_rates_tsvs)],
            error_plot_style = error_plot_style, cell_cutoff_method = cell_cutoff_method, docker = docker, preemptible = preemptible
    }

    call Compare {
        input:
            ready = CheckInputs.ok,
            comparison_name = comparison_name,
            pre_sample_names = pre_sample_names, pre_aligned_length_tsvs = pre_aligned_length_tsvs,
            pre_error_by_position_tsvs = pre_error_by_position_tsvs, pre_read_error_rates_tsvs = pre_read_error_rates_tsvs, pre_cb_counts_tsv_gzs = pre_cb_counts_tsv_gzs,
            pre_summary_jsons = pre_summary_jsons,
            post_sample_names = post_sample_names, post_aligned_length_tsvs = post_aligned_length_tsvs,
            post_dc_tsvs = post_dc_tsvs, post_error_by_position_tsvs = post_error_by_position_tsvs, post_read_error_rates_tsvs = post_read_error_rates_tsvs,
            post_cb_counts_tsv_gzs = post_cb_counts_tsv_gzs, post_summary_jsons = post_summary_jsons,
            error_plot_style = error_plot_style, cell_min_reads = cell_min_reads, cell_cutoff_method = cell_cutoff_method, error_axis_max = error_axis_max, memory_gb = memory_gb, preemptible = preemptible, docker = docker
    }

    output {
        File comparison_summary_tsv = Compare.summary_tsv
        Array[File] comparison_aligned_length_plots = Compare.aligned_length_plots
        Array[File] comparison_dc_plots = Compare.dc_plots
        Array[File] comparison_error_plots = Compare.error_plots
        Array[File] comparison_knee_plots = Compare.knee_plots
        Array[File] comparison_mito_plots = Compare.mito_plots
        Array[File] comparison_combined_pdf = Compare.combined_pdf  # all the figures, one per page
    }
}

task CheckInputs {
    input {
        Int n_pre_names
        Array[Int] n_pre_tables
        Array[Int] n_pre_optional
        Int n_post_names
        Array[Int] n_post_tables
        Array[Int] n_post_optional
        String error_plot_style
        String cell_cutoff_method
        String docker
        Int preemptible
    }

    command <<<
        set -euo pipefail
        err() { echo "ERROR: $1" >&2; exit 1; }
        check() {  # stage, names, required counts (space separated), optional counts
            local stage="$1" names="$2" req="$3" opt="$4" n
            for n in $req; do
                [ "$n" -eq "$names" ] || err "$stage: $names sample names but a table array has $n entries (all required arrays need one entry per name; a sample-name expression that matches no attribute gives an empty list)"
            done
            for n in $opt; do
                [ "$n" -eq 0 ] || [ "$n" -eq "$names" ] || err "$stage: an optional table array has $n entries, expected 0 or $names"
            done
        }
        check pre ~{n_pre_names} "~{sep=" " n_pre_tables}" "~{sep=" " n_pre_optional}"
        check post ~{n_post_names} "~{sep=" " n_post_tables}" "~{sep=" " n_post_optional}"
        [ $(( ~{n_pre_names} + ~{n_post_names} )) -gt 0 ] || err "no samples given: provide the pre_* and/or post_* names and tables"
        case "~{error_plot_style}" in box|violin) ;; *) err "error_plot_style must be box or violin" ;; esac
        case "~{cell_cutoff_method}" in start|steepest) ;; *) err "cell_cutoff_method must be start or steepest" ;; esac
        echo true > ok.txt
    >>>

    output {
        Boolean ok = read_boolean("ok.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 HDD"
        preemptible: preemptible
    }
}

task Compare {
    input {
        Boolean ready  # orders this task after CheckInputs
        String comparison_name
        Array[String] pre_sample_names
        Array[File] pre_aligned_length_tsvs
        Array[File] pre_error_by_position_tsvs
        Array[File] pre_read_error_rates_tsvs
        Array[File] pre_cb_counts_tsv_gzs
        Array[File] pre_summary_jsons
        Array[String] post_sample_names
        Array[File] post_aligned_length_tsvs
        Array[File] post_dc_tsvs
        Array[File] post_error_by_position_tsvs
        Array[File] post_read_error_rates_tsvs
        Array[File] post_cb_counts_tsv_gzs
        Array[File] post_summary_jsons
        String error_plot_style
        Int cell_min_reads
        String cell_cutoff_method
        Float error_axis_max
        Int memory_gb
        Int preemptible
        String docker
    }

    Int disk_gb = ceil(size(pre_cb_counts_tsv_gzs, "GB") * 4 + size(post_cb_counts_tsv_gzs, "GB") * 4 + 10)

    command <<<
        set -euo pipefail
        # inputs validated: ~{ready}
        mkdir -p t out plots
        printf 'sample\tstage\tprefix\n' > manifest.tsv

        pn=(~{sep=" " pre_sample_names}); pa=(~{sep=" " pre_aligned_length_tsvs}); pe=(~{sep=" " pre_error_by_position_tsvs})
        pr=(~{sep=" " pre_read_error_rates_tsvs})
        pc=(~{sep=" " pre_cb_counts_tsv_gzs}); ps=(~{sep=" " pre_summary_jsons})
        for i in "${!pn[@]}"; do
            [ "${#pa[@]}" = "${#pn[@]}" ] && [ "${#pe[@]}" = "${#pn[@]}" ] && [ "${#pc[@]}" = "${#pn[@]}" ] && [ "${#ps[@]}" = "${#pn[@]}" ] \
                || { echo "ERROR: pre arrays differ in length" >&2; exit 1; }
            ln -s "${pa[$i]}" t/pre$i.aligned_length.tsv; ln -s "${pe[$i]}" t/pre$i.error_by_position.tsv
            ln -s "${pc[$i]}" t/pre$i.cb_counts.tsv.gz;   ln -s "${ps[$i]}" t/pre$i.summary.json
            [ "${#pr[@]}" = "${#pn[@]}" ] && ln -s "${pr[$i]}" t/pre$i.read_error_rates.tsv
            printf '%s\tpre\tt/pre%s\n' "${pn[$i]}" "$i" >> manifest.tsv
        done

        qn=(~{sep=" " post_sample_names}); qa=(~{sep=" " post_aligned_length_tsvs}); qd=(~{sep=" " post_dc_tsvs})
        qr=(~{sep=" " post_read_error_rates_tsvs}); qe=(~{sep=" " post_error_by_position_tsvs}); qc=(~{sep=" " post_cb_counts_tsv_gzs}); qs=(~{sep=" " post_summary_jsons})
        for i in "${!qn[@]}"; do
            [ "${#qa[@]}" = "${#qn[@]}" ] && [ "${#qd[@]}" = "${#qn[@]}" ] && [ "${#qe[@]}" = "${#qn[@]}" ] \
                && [ "${#qc[@]}" = "${#qn[@]}" ] && [ "${#qs[@]}" = "${#qn[@]}" ] \
                || { echo "ERROR: post arrays differ in length" >&2; exit 1; }
            ln -s "${qa[$i]}" t/post$i.aligned_length.tsv; ln -s "${qd[$i]}" t/post$i.dc.tsv
            ln -s "${qe[$i]}" t/post$i.error_by_position.tsv; ln -s "${qc[$i]}" t/post$i.cb_counts.tsv.gz
            ln -s "${qs[$i]}" t/post$i.summary.json
            [ "${#qr[@]}" = "${#qn[@]}" ] && ln -s "${qr[$i]}" t/post$i.read_error_rates.tsv
            printf '%s\tpost\tt/post%s\n' "${qn[$i]}" "$i" >> manifest.tsv
        done

        [ "$(wc -l < manifest.tsv)" -gt 1 ] || { echo "ERROR: no samples given" >&2; exit 1; }
        bamqc_compare --manifest manifest.tsv --outdir out --error-style ~{error_plot_style} --cell-min-reads ~{cell_min_reads} --cell-cutoff-method ~{cell_cutoff_method} --error-axis-max ~{error_axis_max}
        for f in out/*; do cp "$f" "plots/~{comparison_name}.$(basename "$f")"; done
    >>>

    output {
        File summary_tsv = "plots/~{comparison_name}.comparison_summary.tsv"
        Array[File] aligned_length_plots = glob("plots/*.comparison_aligned_length_*")
        Array[File] dc_plots = glob("plots/*.comparison_dc.*")
        Array[File] error_plots = glob("plots/*.comparison_error_*")
        Array[File] knee_plots = glob("plots/*.comparison_knee*")
        Array[File] mito_plots = glob("plots/*.comparison_mito*")
        Array[File] combined_pdf = glob("plots/*.comparison_all_plots.pdf")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
