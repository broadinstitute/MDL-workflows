version 1.0

# Plots from the tables written by Bam_Dedup_QC (no BAMs): pre only, post only, or both.
# Rerun it whenever tables are added or changed, e.g. after adding the post-dedup tables to an earlier pre-only run.
# Tables per stage: the four required ones (+ read_error_rates, needed for the per-read error plot; dc for post).
workflow Bam_Dedup_QC_Plot {
    input {
        String sample_name

        File? pre_aligned_length_tsv
        File? pre_error_by_position_tsv
        File? pre_read_error_rates_tsv
        File? pre_cb_counts_tsv_gz
        File? pre_summary_json

        File? post_aligned_length_tsv
        File? post_dc_tsv
        File? post_error_by_position_tsv
        File? post_read_error_rates_tsv
        File? post_cb_counts_tsv_gz
        File? post_summary_json

        Boolean plot_length = true
        Boolean plot_dc = true
        Boolean plot_errors = true
        Boolean plot_knee = true
        String error_plot_style = "box"  # per-read error plot: box or violin

        Int memory_gb = 8
        Int preemptible = 3
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/bam_dedup_qc:0.1.0"
    }

    # A stage is plotted only with its required tables (pre: 4, post: 5 incl. dc); a partial set is an error.
    Int n_pre = length(select_all([pre_aligned_length_tsv, pre_error_by_position_tsv, pre_cb_counts_tsv_gz, pre_summary_json]))
    Int n_post = length(select_all([post_aligned_length_tsv, post_dc_tsv, post_error_by_position_tsv, post_cb_counts_tsv_gz, post_summary_json]))
    Boolean have_pre = n_pre == 4
    Boolean have_post = n_post == 5

    call CheckInputs { input: n_pre = n_pre, n_post = n_post, error_plot_style = error_plot_style, docker = docker, preemptible = preemptible }

    call Plot {
        input:
            ready = CheckInputs.ok,
            sample_name = sample_name,
            pre_aligned_length_tsv = pre_aligned_length_tsv, pre_error_by_position_tsv = pre_error_by_position_tsv,
            pre_read_error_rates_tsv = pre_read_error_rates_tsv, pre_cb_counts_tsv_gz = pre_cb_counts_tsv_gz,
            pre_summary_json = pre_summary_json,
            post_aligned_length_tsv = post_aligned_length_tsv, post_dc_tsv = post_dc_tsv,
            post_error_by_position_tsv = post_error_by_position_tsv, post_read_error_rates_tsv = post_read_error_rates_tsv,
            post_cb_counts_tsv_gz = post_cb_counts_tsv_gz, post_summary_json = post_summary_json,
            have_pre = have_pre, have_post = have_post,
            plot_length = plot_length, plot_dc = plot_dc, plot_errors = plot_errors, plot_knee = plot_knee,
            error_plot_style = error_plot_style, memory_gb = memory_gb, preemptible = preemptible, docker = docker
    }

    output {
        # png + pdf; an empty array = that family was not made (switched off, or its data is missing).
        Array[File] plots_aligned_length = Plot.plots_aligned_length
        Array[File] plots_dc = Plot.plots_dc
        Array[File] plots_error_rates = Plot.plots_error_rates
        Array[File] plots_knee = Plot.plots_knee
        File plots_stats_json = Plot.plots_stats_json
    }
}

task CheckInputs {
    input {
        Int n_pre
        Int n_post
        String error_plot_style
        String docker
        Int preemptible
    }

    command <<<
        set -euo pipefail
        err() { echo "ERROR: $1" >&2; exit 1; }
        [ ~{n_pre} -eq 0 ] || [ ~{n_pre} -eq 4 ] || err "~{n_pre} of the 4 required pre_* tables given (aligned_length, error_by_position, cb_counts, summary)"
        [ ~{n_post} -eq 0 ] || [ ~{n_post} -eq 5 ] || err "~{n_post} of the 5 required post_* tables given (aligned_length, dc, error_by_position, cb_counts, summary)"
        [ ~{n_pre} -gt 0 ] || [ ~{n_post} -gt 0 ] || err "no tables given: provide the pre_* and/or the post_* tables"
        case "~{error_plot_style}" in box|violin) ;; *) err "error_plot_style must be box or violin" ;; esac
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

task Plot {
    input {
        Boolean ready  # orders this task after CheckInputs
        String sample_name
        File? pre_aligned_length_tsv
        File? pre_error_by_position_tsv
        File? pre_read_error_rates_tsv
        File? pre_cb_counts_tsv_gz
        File? pre_summary_json
        File? post_aligned_length_tsv
        File? post_dc_tsv
        File? post_error_by_position_tsv
        File? post_read_error_rates_tsv
        File? post_cb_counts_tsv_gz
        File? post_summary_json
        Boolean have_pre
        Boolean have_post
        Boolean plot_length
        Boolean plot_dc
        Boolean plot_errors
        Boolean plot_knee
        String error_plot_style
        Int memory_gb
        Int preemptible
        String docker
    }

    String families = (if plot_length then "length," else "") + (if plot_dc then "dc," else "") +
                      (if plot_errors then "errors," else "") + (if plot_knee then "knee," else "")
    Int disk_gb = ceil(size(pre_cb_counts_tsv_gz, "GB") * 4 + size(post_cb_counts_tsv_gz, "GB") * 4 + 10)

    command <<<
        set -euo pipefail
        # inputs validated: ~{ready}
        mkdir -p t raw plots
        args=""
        if [ "~{have_pre}" = "true" ]; then
            ln -s ~{pre_aligned_length_tsv} t/pre.aligned_length.tsv
            ln -s ~{pre_error_by_position_tsv} t/pre.error_by_position.tsv
            ln -s ~{pre_cb_counts_tsv_gz} t/pre.cb_counts.tsv.gz
            ln -s ~{pre_summary_json} t/pre.summary.json
            if [ -n "~{pre_read_error_rates_tsv}" ]; then ln -s ~{pre_read_error_rates_tsv} t/pre.read_error_rates.tsv; fi
            args="$args --pre t/pre"
        fi
        if [ "~{have_post}" = "true" ]; then
            ln -s ~{post_aligned_length_tsv} t/post.aligned_length.tsv
            ln -s ~{post_dc_tsv} t/post.dc.tsv
            ln -s ~{post_error_by_position_tsv} t/post.error_by_position.tsv
            ln -s ~{post_cb_counts_tsv_gz} t/post.cb_counts.tsv.gz
            ln -s ~{post_summary_json} t/post.summary.json
            if [ -n "~{post_read_error_rates_tsv}" ]; then ln -s ~{post_read_error_rates_tsv} t/post.read_error_rates.tsv; fi
            args="$args --post t/post"
        fi
        bamqc_plot $args --plots "~{families}" --error-style ~{error_plot_style} --outdir raw
        for f in raw/*; do cp "$f" "plots/~{sample_name}.$(basename "$f")"; done
    >>>

    output {
        Array[File] plots_aligned_length = glob("plots/*.aligned_length_*")
        Array[File] plots_dc = glob("plots/*.dc_distribution.*")
        Array[File] plots_error_rates = glob("plots/*.error_rates_*")
        Array[File] plots_knee = glob("plots/*.cell_knee.*")
        File plots_stats_json = "plots/~{sample_name}.qc_stats.json"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible
    }
}
