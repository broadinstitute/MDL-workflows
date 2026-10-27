version 1.1

task BC_Correct_Pass1 {
    input {
        File input_fastq
        String sample_name
        String whitelist
        String bc_range
        String umi_range
        Int trim_len
        String end_adapter
        String start_adapter
        String endedness
        Float min_post
        Int min_trimmed_len
        Boolean qual_tags = true
        Boolean revcomp_trimmed = false
        Boolean rescue_start_adapter = false
        String docker

        String poly_logic = ""
        String poly_start = ""
        String poly_window = ""
        String poly_t_min = ""
        String poly_a_min = ""
        String poly_err = ""

        Int? disk_size_gb
    }

    Boolean poly_enabled = (poly_logic != "") || (poly_start != "") || (poly_window != "") || (poly_t_min != "") || (poly_a_min != "") || (poly_err != "")
    String resolved_poly_logic = if (poly_enabled) then (if (poly_logic != "") then poly_logic else "cutadapt") else ""
    String poly_start_arg = if (poly_start != "") then ("--poly-start " + poly_start) else ""
    String poly_window_arg = if (poly_window != "") then ("--poly-window " + poly_window) else ""
    String poly_t_min_arg = if (poly_t_min != "") then ("--polyT-min " + poly_t_min) else ""
    String poly_a_min_arg = if (poly_a_min != "") then ("--polyA-min " + poly_a_min) else ""
    String poly_err_arg = if (poly_err != "") then ("--poly-err " + poly_err) else ""
    String poly_logic_arg = if (poly_enabled) then ("--poly-logic " + resolved_poly_logic) else ""
    String qual_tags_arg = if (qual_tags) then "--qual-tags" else ""
    String rescue_arg = if (rescue_start_adapter) then "--rescue-start-adapter" else ""
    String revcomp_arg = if (revcomp_trimmed) then "--revcomp-trimmed" else ""

    String trimmed_prefix = sample_name
    String counts_table_name = trimmed_prefix + ".whitelist.counts.txt.gz"
    String matched_fastq_path = trimmed_prefix + ".matched.fastq.gz"
    String unmatched_fastq_path = trimmed_prefix + ".unmatched.fastq.gz"

    Int diskGB = select_first([disk_size_gb, ceil(size(input_fastq, "GB") * 3 + 20)])

    command <<<
        set -euo pipefail

        if [[ -z "~{whitelist}" ]]; then
          echo "Error: whitelist is required. Provide library_config (resolved whitelist) or whitelist input." >&2
          exit 1
        fi

        /opt/pipeline/fastq_preprocessing_pipeline/trim_split_nocap_extracttags.sh \
          --stdout-found \
          ~{rescue_arg} \
          -i ~{input_fastq} \
          -o ~{trimmed_prefix} \
          -t 2 \
          -a ~{end_adapter} \
          -g ~{start_adapter} \
          -e ~{endedness} \
        | bc_correct_split \
          --fastq4 \
          --whitelist ~{whitelist} \
          --outdir . \
          --threads 2 \
          --bc ~{bc_range} \
          --umi ~{umi_range} \
          --trim ~{trim_len} \
          --min-post ~{min_post} \
          --min-trimmed-len ~{min_trimmed_len} \
          --pass 1 \
          --counts-table ~{counts_table_name} \
          --output-prefix ~{trimmed_prefix} \
          ~{qual_tags_arg} \
          ~{revcomp_arg} \
          ~{poly_start_arg} \
          ~{poly_window_arg} \
          ~{poly_t_min_arg} \
          ~{poly_a_min_arg} \
          ~{poly_err_arg} \
          ~{poly_logic_arg} \
          -
        tail -n +2 "~{trimmed_prefix}_trim_stats.txt" | tr '\t\n' '=;' > "~{trimmed_prefix}_trim_stats.txt.oneline"
        tail -n +2 "~{trimmed_prefix}.pass1.readstats.txt" | tr '\t\n' '=;' > "~{trimmed_prefix}.pass1.readstats.txt.oneline"
    >>>

    output {
        File matched_fastq = "~{matched_fastq_path}"
        File unmatched_fastq = "~{unmatched_fastq_path}"
        File counts_table = "~{counts_table_name}"
        File trim_stats = "~{trimmed_prefix}_trim_stats.txt"
        File read_stats = "~{trimmed_prefix}.pass1.readstats.txt"
        # One-line "key=value;..." copies so the accounting task needs no per-file localization.
        String trim_stats_line = read_string("~{trimmed_prefix}_trim_stats.txt.oneline")
        String read_stats_line = read_string("~{trimmed_prefix}.pass1.readstats.txt.oneline")
    }

    # cpu/memory intentionally omitted: GCP Batch has a known bug where specifying both
    # predefinedMachineType and an explicit compute_resource (cpu_milli/memory_mib) can spuriously
    # reject an otherwise-valid combination ("machine_type ... cannot satisfy compute_resource
    # ..."), even when the values exactly match the machine type's real spec. Workaround: let
    # predefinedMachineType alone determine the shape (c3d-highcpu-4 = 4 vCPU / 8 GB).
    runtime {
        docker: docker
        predefinedMachineType: "c3d-highcpu-4"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}



task Merge_Whitelist_Counts {
    input {
        Array[File] counts_tables
        String docker
        String output_path = "merged.whitelist.counts.txt.gz"
        Boolean gzip_output = true
        Int gzip_level = 1
        Int? disk_size_gb
    }

    String effective_script = "/opt/pipeline/fastq_preprocessing_pipeline/sum_whitelist_counts.py"
    Int diskGB = select_first([disk_size_gb, ceil(size(counts_tables, "GB") * 3 + 20)])
    String gzip_output_arg = if (gzip_output) then "--gzip-output --gzip-level " + gzip_level else ""

    command <<<
        set -euo pipefail
        all_counts=(~{sep=" " counts_tables})
        python3 ~{effective_script} ~{gzip_output_arg} --output ~{output_path} ${all_counts[@]}
    >>>

    output {
        File merged_counts = "~{output_path}"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "4 GB"
        predefinedMachineType: "n2d-custom-2-4096"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Resolve_BC_Params_From_Config {
    input {
        String library_config
        String docker
    }

    command <<<
        set -euo pipefail

        value="~{library_config}"
        if [[ "$value" == gs://* ]]; then
          gsutil cp "$value" ./custom_library_config.yaml
          config_path="./custom_library_config.yaml"
        elif [[ "$value" == *.yaml ]] || [[ "$value" == *.yml ]]; then
          config_path="$value"
        else
          preset_base="${value##*/}"
          case "$preset_base" in
            "10x_3p_v3"|"10x_3p_v3.1")
              preset_base="10x_3p_v3-3.1"
              ;;
          esac
          config_path="/opt/pipeline/fastq_preprocessing_pipeline/configs/${preset_base}.yaml"
        fi
        if [ ! -f "$config_path" ]; then
          echo "Error: config file not found at $config_path" >&2
          exit 1
        fi
        export CONFIG_PATH="$config_path"

        python3 - <<'PY'
        import importlib.util
        import os

        spec = importlib.util.spec_from_file_location(
            "pipeline_wrapper", "/opt/pipeline/fastq_preprocessing_pipeline/pipeline_wrapper.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        config_path = os.environ["CONFIG_PATH"]
        config = module.load_config(config_path)
        params = module.build_from_config(config)

        def write(path, value, default="50"):
            with open(path, "w", encoding="utf-8") as f:
                if value is None:
                    f.write(str(default))
                elif isinstance(value, bool):
                    f.write(str(value).lower())
                else:
                    f.write(str(value))

        write("resolved_whitelist.txt", params["whitelist"])
        write("resolved_endedness.txt", params["endedness"])
        write("resolved_start_adapter.txt", params["start_adapter"])
        write("resolved_end_adapter.txt", params["end_adapter"])
        write("resolved_bc_range.txt", params["bc_range"])
        write("resolved_umi_range.txt", params["umi_range"])
        write("resolved_trim_len.txt", params["trim_len"])
        write("resolved_min_trimmed_len.txt", params["min_trimmed_len"], default="50")
        write("resolved_revcomp_trimmed.txt", params.get("cdna_reverse_complement"), default="false")
        poly = params.get("poly")
        if not isinstance(poly, dict):
            poly = {}
        write("resolved_poly_logic.txt", poly.get("logic"), default="")
        write("resolved_poly_start.txt", poly.get("start"), default="")
        write("resolved_poly_window.txt", poly.get("window"), default="")
        write("resolved_poly_t_min.txt", poly.get("min_t"), default="")
        write("resolved_poly_a_min.txt", poly.get("min_a"), default="")
        write("resolved_poly_err.txt", poly.get("err"), default="")
        PY
    >>>

    output {
        String resolved_whitelist = read_string("resolved_whitelist.txt")
        String resolved_endedness = read_string("resolved_endedness.txt")
        String resolved_start_adapter = read_string("resolved_start_adapter.txt")
        String resolved_end_adapter = read_string("resolved_end_adapter.txt")
        String resolved_bc_range = read_string("resolved_bc_range.txt")
        String resolved_umi_range = read_string("resolved_umi_range.txt")
        Int resolved_trim_len = read_int("resolved_trim_len.txt")
        Int resolved_min_trimmed_len = read_int("resolved_min_trimmed_len.txt")
        Boolean resolved_revcomp_trimmed = read_boolean("resolved_revcomp_trimmed.txt")
        String resolved_poly_logic = read_string("resolved_poly_logic.txt")
        String resolved_poly_start = read_string("resolved_poly_start.txt")
        String resolved_poly_window = read_string("resolved_poly_window.txt")
        String resolved_poly_t_min = read_string("resolved_poly_t_min.txt")
        String resolved_poly_a_min = read_string("resolved_poly_a_min.txt")
        String resolved_poly_err = read_string("resolved_poly_err.txt")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        predefinedMachineType: "n1-custom-1-1024"
        preemptible: 3
    }
}


task BC_Correct_Pass2 {
    input {
        File unmatched_fastq
        File counts_table
        String sample_name
        String whitelist
        String bc_range
        String umi_range
        Int trim_len
        Float min_post
        Int min_trimmed_len
        Boolean qual_tags = true
        Boolean revcomp_trimmed = false
        Boolean generate_discarded = false
        String docker

        String poly_logic = ""
        String poly_start = ""
        String poly_window = ""
        String poly_t_min = ""
        String poly_a_min = ""
        String poly_err = ""

        Int? disk_size_gb
    }

    Boolean poly_enabled = (poly_logic != "") || (poly_start != "") || (poly_window != "") || (poly_t_min != "") || (poly_a_min != "") || (poly_err != "")
    String resolved_poly_logic = if (poly_enabled) then (if (poly_logic != "") then poly_logic else "cutadapt") else ""
    String poly_start_arg = if (poly_start != "") then ("--poly-start " + poly_start) else ""
    String poly_window_arg = if (poly_window != "") then ("--poly-window " + poly_window) else ""
    String poly_t_min_arg = if (poly_t_min != "") then ("--polyT-min " + poly_t_min) else ""
    String poly_a_min_arg = if (poly_a_min != "") then ("--polyA-min " + poly_a_min) else ""
    String poly_err_arg = if (poly_err != "") then ("--poly-err " + poly_err) else ""
    String poly_logic_arg = if (poly_enabled) then ("--poly-logic " + resolved_poly_logic) else ""
    String qual_tags_arg = if (qual_tags) then "--qual-tags" else ""
    String revcomp_arg = if (revcomp_trimmed) then "--revcomp-trimmed" else ""
    String discard_arg = if (generate_discarded) then "--discarded" else "--no-discarded"

    String output_prefix = sample_name
    String corrected_fastq_path = output_prefix + ".corrected.fastq.gz"

    Int diskGB = select_first([disk_size_gb, ceil(size(unmatched_fastq, "GB") * 3 + 20)])

    command <<<
        set -euo pipefail

        if [[ -z "~{whitelist}" ]]; then
          echo "Error: whitelist is required. Provide library_config (resolved whitelist) or whitelist input." >&2
          exit 1
        fi

        bc_correct_split \
          --fastq4 \
          --whitelist ~{whitelist} \
          --outdir . \
          --threads 2 \
          --bc ~{bc_range} \
          --umi ~{umi_range} \
          --trim ~{trim_len} \
          --min-post ~{min_post} \
          --min-trimmed-len ~{min_trimmed_len} \
          --pass 2 \
          --counts-table ~{counts_table} \
          --output-prefix ~{output_prefix} \
          ~{discard_arg} \
          ~{qual_tags_arg} \
          ~{revcomp_arg} \
          ~{poly_start_arg} \
          ~{poly_window_arg} \
          ~{poly_t_min_arg} \
          ~{poly_a_min_arg} \
          ~{poly_err_arg} \
          ~{poly_logic_arg} \
          ~{unmatched_fastq}
        tail -n +2 "~{output_prefix}.pass2.readstats.txt" | tr '\t\n' '=;' > "~{output_prefix}.pass2.readstats.txt.oneline"
    >>>

    output {
        File corrected_fastq = "~{corrected_fastq_path}"
        Array[File] discarded_fastqs = glob("~{output_prefix}*.discarded*.fastq.gz")
        File read_stats = "~{output_prefix}.pass2.readstats.txt"
        String read_stats_line = read_string("~{output_prefix}.pass2.readstats.txt.oneline")
    }

    # cpu/memory intentionally omitted: GCP Batch has a known bug where specifying both
    # predefinedMachineType and an explicit compute_resource (cpu_milli/memory_mib) can spuriously
    # reject an otherwise-valid combination ("machine_type ... cannot satisfy compute_resource
    # ..."), even when the values exactly match the machine type's real spec. Workaround: let
    # predefinedMachineType alone determine the shape (c3d-highcpu-4 = 4 vCPU / 8 GB).
    runtime {
        docker: docker
        predefinedMachineType: "c3d-highcpu-4"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Summarize_Read_Accounting {
    # Per-shard and pool-total read accounting from the pass1/pass2 stats files, so the reads that
    # are lost at each filter (trimming, too short, uncorrectable) and the reads rescued by error
    # correction are tracked. exact_match_written is the exact-whitelist CB count from pass 1 (the
    # same reads that feed the merged counts table).
    input {
        Array[String] sample_names
        Array[String] trim_stats
        Array[String] pass1_read_stats
        Array[String] pass2_read_stats
        String docker
    }

    command <<<
        set -euo pipefail

        python3 - <<'PY'
        names = ["~{sep='", "' sample_names}"]
        trim = open("~{write_lines(trim_stats)}").read().splitlines()
        p1 = open("~{write_lines(pass1_read_stats)}").read().splitlines()
        p2 = open("~{write_lines(pass2_read_stats)}").read().splitlines()

        def load(line):
            return {k: int(v) for k, v in (kv.split("=") for kv in line.strip(";").split(";"))}

        cols = [
            ("raw_reads", lambda t, a, b: t["stage1_input_reads"]),
            ("trim_dropped_stage1", lambda t, a, b: t["stage1_dropped_reads"]),
            ("trim_dropped_stage2", lambda t, a, b: t["stage2_dropped_reads"]),
            ("reads_into_bc_correct", lambda t, a, b: a["total_reads"]),
            ("exact_match_written", lambda t, a, b: a["matched_written"]),
            ("exact_match_too_short_dropped", lambda t, a, b: a["matched_too_short_dropped"]),
            ("unmatched_to_pass2", lambda t, a, b: a["unmatched_written"]),
            ("rescued_by_correction", lambda t, a, b: b["corrected_written"]),
            ("corrected_too_short_dropped", lambda t, a, b: b["dropped_too_short"]),
            ("uncorrectable_dropped", lambda t, a, b: b["dropped_uncorrectable"]),
            ("out_of_range_dropped", lambda t, a, b: b["dropped_out_of_range"]),
            ("written_to_discarded", lambda t, a, b: b["written_to_discarded"]),
            ("total_kept", lambda t, a, b: a["matched_written"] + b["corrected_written"]),
            # start-adapter rescue (0 when the run did not use it); these reads are part of the outcomes above, not extra
            ("start_rescued_substitution", lambda t, a, b: t.get("rescued_substitution_reads", 0)),
            ("start_rescued_indel", lambda t, a, b: t.get("rescued_indel_reads", 0)),
            # The same outcomes counted for the rescued reads only (XR tag), by bc_correct_split itself at each decision point
            # (0 when the run did not use the rescue or the tool predates these counters).
            ("xr_exact_match_written", lambda t, a, b: a.get("matched_written_xr", 0)),
            ("xr_exact_match_too_short_dropped", lambda t, a, b: a.get("matched_too_short_dropped_xr", 0)),
            ("xr_unmatched_to_pass2", lambda t, a, b: a.get("unmatched_written_xr", 0)),
            ("xr_rescued_by_correction", lambda t, a, b: b.get("corrected_written_xr", 0)),
            ("xr_corrected_too_short_dropped", lambda t, a, b: b.get("dropped_too_short_xr", 0)),
            ("xr_uncorrectable_dropped", lambda t, a, b: b.get("dropped_uncorrectable_xr", 0)),
            ("xr_out_of_range_dropped", lambda t, a, b: b.get("dropped_out_of_range_xr", 0)),
            ("xr_written_to_discarded", lambda t, a, b: b.get("written_to_discarded_xr", 0)),
            ("xr_total_kept", lambda t, a, b: a.get("matched_written_xr", 0) + b.get("corrected_written_xr", 0)),
        ]
        header = ["sample_name"] + [c for c, _ in cols]

        totals = [0] * len(cols)
        with open("read_accounting_per_shard.tsv", "w") as out:
            out.write("\t".join(header) + "\n")
            for name, tp, ap, bp in zip(names, trim, p1, p2):
                t, a, b = load(tp), load(ap), load(bp)
                vals = [fn(t, a, b) for _, fn in cols]
                totals = [x + y for x, y in zip(totals, vals)]
                out.write("\t".join([name] + [str(v) for v in vals]) + "\n")

        with open("read_accounting_total.tsv", "w") as out:
            out.write("metric\tcount\tfraction_of_raw_reads\n")
            for (c, _), v in zip(cols, totals):
                frac = v / totals[0] if totals[0] else 0
                out.write(f"{c}\t{v}\t{frac:.6f}\n")

        # What the run would have looked like without the start-adapter rescue: every rescued read (counted at trimming) is back in
        # trim_dropped_stage2 and out of the later outcomes (its fate there is the xr_* counters).
        tot = {c: v for (c, _), v in zip(cols, totals)}
        rescued = tot["start_rescued_substitution"] + tot["start_rescued_indel"]
        no_rescue = [
            ("raw_reads", tot["raw_reads"]),
            ("trim_dropped_stage1", tot["trim_dropped_stage1"]),
            ("trim_dropped_stage2", tot["trim_dropped_stage2"] + rescued),
            ("reads_into_bc_correct", tot["reads_into_bc_correct"] - rescued),
            ("exact_match_written", tot["exact_match_written"] - tot["xr_exact_match_written"]),
            ("exact_match_too_short_dropped", tot["exact_match_too_short_dropped"] - tot["xr_exact_match_too_short_dropped"]),
            ("unmatched_to_pass2", tot["unmatched_to_pass2"] - tot["xr_unmatched_to_pass2"]),
            ("rescued_by_correction", tot["rescued_by_correction"] - tot["xr_rescued_by_correction"]),
            ("corrected_too_short_dropped", tot["corrected_too_short_dropped"] - tot["xr_corrected_too_short_dropped"]),
            ("uncorrectable_dropped", tot["uncorrectable_dropped"] - tot["xr_uncorrectable_dropped"]),
            ("out_of_range_dropped", tot["out_of_range_dropped"] - tot["xr_out_of_range_dropped"]),
            ("written_to_discarded", tot["written_to_discarded"] - tot["xr_written_to_discarded"]),
            ("total_kept", tot["total_kept"] - tot["xr_total_kept"]),
        ]
        with open("read_accounting_total_without_rescue.tsv", "w") as out:
            out.write("metric\tcount\tfraction_of_raw_reads\n")
            for c, v in no_rescue:
                out.write(f"{c}\t{v}\t{(v / tot['raw_reads'] if tot['raw_reads'] else 0):.6f}\n")
        PY
    >>>

    output {
        File per_shard = "read_accounting_per_shard.tsv"
        File total = "read_accounting_total.tsv"
        File total_without_rescue = "read_accounting_total_without_rescue.tsv"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "1 GB"
        predefinedMachineType: "n1-custom-1-1024"
        preemptible: 3
    }
}


task Plot_Read_Accounting {
    # Plots of the Summarize_Read_Accounting tables (numpy + matplotlib only): page 1 is where the raw
    # reads of the whole set end up, page 2 the same per shard (to spot shards that behave differently).
    # "unaccounted" is raw_reads minus every listed outcome; it should be 0 and is drawn if it is not.
    input {
        File per_shard
        File total
        String plot_label
        String docker
    }

    command <<<
        set -euo pipefail

        python3 - <<'PY'
        import csv
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        label = "~{plot_label}"
        KEPT, LOST, ODD = "#2c7fb8", "#d95f0e", "#7f7f7f"

        with open("~{total}") as f:
            tot = {r["metric"]: int(r["count"]) for r in csv.DictReader(f, delimiter="\t")}
        with open("~{per_shard}") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        names = [r["sample_name"] for r in rows]
        col = lambda c: np.array([int(r[c]) for r in rows], dtype=float)

        raw = tot["raw_reads"]
        trim_drop = tot["trim_dropped_stage1"] + tot["trim_dropped_stage2"]
        outcomes = [
            ("trimming: dropped (adapters / length)", trim_drop, LOST),
            ("exact barcode: too short after trimming", tot["exact_match_too_short_dropped"], LOST),
            ("exact barcode: kept", tot["exact_match_written"], KEPT),
            ("corrected barcode: kept (rescued)", tot["rescued_by_correction"], KEPT),
            ("corrected barcode: too short", tot["corrected_too_short_dropped"], LOST),
            ("barcode uncorrectable", tot["uncorrectable_dropped"], LOST),
            ("barcode out of range", tot["out_of_range_dropped"], LOST),
        ]
        unaccounted = raw - sum(c for _, c, _ in outcomes)
        if unaccounted != 0:
            outcomes.append(("UNACCOUNTED (check the stats files)", unaccounted, ODD))
        kept = tot["total_kept"]
        den = max(raw, 1)
        resc_sub, resc_ind = tot.get("start_rescued_substitution", 0), tot.get("start_rescued_indel", 0)
        title = f"{label}: " if label else ""

        with PdfPages("read_accounting_plots.pdf") as pdf:
            fig, ax = plt.subplots(figsize=(11, 5.5))
            y = np.arange(len(outcomes))[::-1]
            ax.barh(y, [100 * c / den for _, c, _ in outcomes], color=[k for _, _, k in outcomes])
            ax.set_yticks(y)
            ax.set_yticklabels([n for n, _, _ in outcomes])
            top = max([100 * c / den for _, c, _ in outcomes] + [1])
            for yi, (_, c, _) in zip(y, outcomes):
                ax.text(100 * c / den + top * 0.01, yi, f"{c:,}  ({100 * c / den:.2f}%)", va="center", fontsize=9)
            ax.set_xlim(0, top * 1.35)
            ax.set_xlabel("% of raw reads")
            ax.set_title(f"{title}{100 * kept / den:.1f}% of {raw:,} raw reads kept ({len(rows)} shards)")
            ax.spines[["top", "right"]].set_visible(False)
            if resc_sub + resc_ind > 0:
                ax.text(0.99, 0.02, f"start-adapter rescue (included above): {resc_sub:,} substitution + {resc_ind:,} indel reads = {100 * (resc_sub + resc_ind) / den:.2f}% of raw reads",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#444444")
            fig.tight_layout()
            fig.savefig("read_accounting_summary.png", dpi=150)
            pdf.savefig(fig)
            plt.close(fig)

            raw_s = col("raw_reads")
            d = np.maximum(raw_s, 1)
            parts = [
                ("exact barcode: kept", col("exact_match_written"), KEPT),
                ("corrected barcode: kept", col("rescued_by_correction"), "#7fcdbb"),
                ("trimming: dropped", col("trim_dropped_stage1") + col("trim_dropped_stage2"), "#fdae61"),
                ("too short after trimming", col("exact_match_too_short_dropped") + col("corrected_too_short_dropped"), "#d7191c"),
                ("uncorrectable / out of range", col("uncorrectable_dropped") + col("out_of_range_dropped"), "#7b3294"),
            ]
            kept_frac = col("total_kept") / d
            x = np.arange(len(rows))
            fig, axs = plt.subplots(3, 1, figsize=(11, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1.5, 1.5]})
            bottom = np.zeros(len(rows))
            for n, v, c in parts:
                axs[0].bar(x, v / d, bottom=bottom, width=1.0, color=c, label=n, linewidth=0)
                bottom += v / d
            axs[0].set_ylabel("fraction of the shard's raw reads")
            axs[0].set_ylim(0, 1)
            axs[0].legend(loc="lower right", fontsize=8, ncol=2, framealpha=0.9)
            axs[0].set_title(f"{title}per-shard outcome of the raw reads")
            axs[1].bar(x, raw_s / 1e6, width=1.0, color="#636363", linewidth=0)
            axs[1].set_ylabel("raw reads (M)")
            axs[2].plot(x, kept_frac, ".", ms=3, color=KEPT)
            med = float(np.median(kept_frac))
            mad = float(np.median(np.abs(kept_frac - med)))
            tol = max(0.05, 5 * mad)
            axs[2].axhline(med, color="k", lw=0.8)
            axs[2].axhspan(med - tol, med + tol, color="k", alpha=0.08)
            out = [i for i in range(len(rows)) if abs(kept_frac[i] - med) > tol]
            axs[2].plot(x[out], kept_frac[out], "o", ms=5, mfc="none", color=LOST)
            axs[2].set_ylabel("kept fraction")
            axs[2].set_xlabel("shard (input order)")
            axs[2].set_xlim(-0.5, len(rows) - 0.5)
            axs[2].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            for a in axs:
                a.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            lines = [f"{names[i]}   kept {kept_frac[i]:.3f}   raw reads {int(raw_s[i]):,}" for i in out[:25]]
            if len(out) > 25:
                lines.append(f"... and {len(out) - 25} more")
            print(f"{len(out)} shard(s) with kept fraction more than {tol:.3f} from the median {med:.3f}")
            if out:
                fig = plt.figure(figsize=(11, 4))
                fig.text(0.02, 0.95, f"{title}shards whose kept fraction differs from the median ({med:.3f}) by more than {tol:.3f}", va="top", fontsize=10)
                fig.text(0.02, 0.85, "\n".join(lines), va="top", fontsize=8, family="monospace")
                pdf.savefig(fig)
                plt.close(fig)
        PY
    >>>

    output {
        File pdf = "read_accounting_plots.pdf"
        File summary_png = "read_accounting_summary.png"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "2 GB"
        predefinedMachineType: "n1-custom-1-2048"
        preemptible: 3
    }
}


workflow BC_Barcode_Extract_And_Correct_Array {
    meta {
        description: "Extract and correct cell barcodes/UMIs in two passes over arrays of FASTQs."
    }

    input {
        Array[File]+ raw_fastqs
        Array[String] sample_names
        File? whitelist
        Float min_post = 0.975
        String library_config

        Boolean qual_tags = true
        Boolean generate_discarded = false
        Boolean output_unmatched = false

        # Rescue of reads that fail the start-adapter filter because of a few edits in the left (TSO/RT) or right (SLS) adapter
        # (anchored, CB+UMI window must keep its exact length). Off by default; needs an image with stage2_rescue; validated on 5' and 3' 10x presets only.
        Boolean rescue_start_adapter = false

        File? whitelist_override

        String docker

        # Plots of the read accounting: needs python3 with numpy + matplotlib (the barcode image has neither).
        String plot_docker = "us-central1-docker.pkg.dev/methods-dev-lab/misc-utilities/bam_dedup_qc:0.1.0"
        String plot_label = ""
    }

    # The counts merge below is a fixed 5-round tree (merge_chunk_size per chunk); the check here guards its capacity.
    Int merge_chunk_size = 10
    Int merge_tree_capacity = merge_chunk_size * merge_chunk_size * merge_chunk_size * merge_chunk_size * merge_chunk_size

    # Fails the workflow (no task) when raw_fastqs is empty or longer than the merge tree can collapse; select_first of an unset optional is
    # only evaluated in that branch. (Cromwell also rejects an empty raw_fastqs earlier because of the Array[File]+ type.)
    File? no_file = None
    Array[File] error_if_fastqs_exceed_merge_capacity = if (length(raw_fastqs) < 1 || length(raw_fastqs) > merge_tree_capacity) then [select_first([no_file])] else raw_fastqs

    Array[Int] lane_indexes = range(length(error_if_fastqs_exceed_merge_capacity))

    call Resolve_BC_Params_From_Config as resolved_config {
        input:
            library_config = library_config,
            docker = docker
    }

    String resolved_whitelist_path = "/opt/pipeline/" + resolved_config.resolved_whitelist
    String effective_whitelist = select_first([whitelist_override, whitelist, resolved_whitelist_path])
    String effective_endedness = resolved_config.resolved_endedness
    String effective_start_adapter = resolved_config.resolved_start_adapter
    String effective_end_adapter = resolved_config.resolved_end_adapter
    String effective_bc_range = resolved_config.resolved_bc_range
    String effective_umi_range = resolved_config.resolved_umi_range
    Int effective_trim_len = resolved_config.resolved_trim_len
    Int effective_min_trimmed_len = resolved_config.resolved_min_trimmed_len

    scatter (i in lane_indexes) {
        call BC_Correct_Pass1 as pass1 {
            input:
                input_fastq = error_if_fastqs_exceed_merge_capacity[i],
                sample_name = sample_names[i],
                whitelist = effective_whitelist,
                bc_range = effective_bc_range,
                umi_range = effective_umi_range,
                trim_len = effective_trim_len,
                end_adapter = effective_end_adapter,
                start_adapter = effective_start_adapter,
                endedness = effective_endedness,
                min_post = min_post,
                min_trimmed_len = effective_min_trimmed_len,
                qual_tags = qual_tags,
                revcomp_trimmed = resolved_config.resolved_revcomp_trimmed,
                poly_logic = resolved_config.resolved_poly_logic,
                poly_start = resolved_config.resolved_poly_start,
                poly_window = resolved_config.resolved_poly_window,
                poly_t_min = resolved_config.resolved_poly_t_min,
                poly_a_min = resolved_config.resolved_poly_a_min,
                poly_err = resolved_config.resolved_poly_err,
                rescue_start_adapter = rescue_start_adapter,
                docker = docker
        }
    }

    Int round1_total_chunks = length(pass1.counts_table)
    Int round1_num_chunks = (round1_total_chunks + merge_chunk_size - 1) / merge_chunk_size
    Array[Int] round1_indexes = range(round1_num_chunks)

    scatter (r1_idx in round1_indexes) {
        Int r1_file_idx_0 = r1_idx * merge_chunk_size + 0
        Int r1_file_idx_1 = r1_idx * merge_chunk_size + 1
        Int r1_file_idx_2 = r1_idx * merge_chunk_size + 2
        Int r1_file_idx_3 = r1_idx * merge_chunk_size + 3
        Int r1_file_idx_4 = r1_idx * merge_chunk_size + 4
        Int r1_file_idx_5 = r1_idx * merge_chunk_size + 5
        Int r1_file_idx_6 = r1_idx * merge_chunk_size + 6
        Int r1_file_idx_7 = r1_idx * merge_chunk_size + 7
        Int r1_file_idx_8 = r1_idx * merge_chunk_size + 8
        Int r1_file_idx_9 = r1_idx * merge_chunk_size + 9
        Array[File] r1_chunk = select_all([
            if (r1_file_idx_0 < round1_total_chunks) then pass1.counts_table[r1_file_idx_0] else None,
            if (r1_file_idx_1 < round1_total_chunks) then pass1.counts_table[r1_file_idx_1] else None,
            if (r1_file_idx_2 < round1_total_chunks) then pass1.counts_table[r1_file_idx_2] else None,
            if (r1_file_idx_3 < round1_total_chunks) then pass1.counts_table[r1_file_idx_3] else None,
            if (r1_file_idx_4 < round1_total_chunks) then pass1.counts_table[r1_file_idx_4] else None,
            if (r1_file_idx_5 < round1_total_chunks) then pass1.counts_table[r1_file_idx_5] else None,
            if (r1_file_idx_6 < round1_total_chunks) then pass1.counts_table[r1_file_idx_6] else None,
            if (r1_file_idx_7 < round1_total_chunks) then pass1.counts_table[r1_file_idx_7] else None,
            if (r1_file_idx_8 < round1_total_chunks) then pass1.counts_table[r1_file_idx_8] else None,
            if (r1_file_idx_9 < round1_total_chunks) then pass1.counts_table[r1_file_idx_9] else None
        ])

        call Merge_Whitelist_Counts as r1_merge {
            input:
                counts_tables = r1_chunk,
                output_path = "merged.whitelist.counts.r1." + r1_idx + ".txt.gz",
                gzip_output = true,
                gzip_level = 1,
                docker = docker
        }
    }

    Array[File] round1_counts = r1_merge.merged_counts

    Int round2_total_chunks = length(round1_counts)
    # round 2 only runs when the previous round left more than one file
    if (round2_total_chunks > 1) {
        Int round2_num_chunks = (round2_total_chunks + merge_chunk_size - 1) / merge_chunk_size
        Array[Int] round2_indexes = range(round2_num_chunks)

        scatter (r2_idx in round2_indexes) {
            Int r2_file_idx_0 = r2_idx * merge_chunk_size + 0
            Int r2_file_idx_1 = r2_idx * merge_chunk_size + 1
            Int r2_file_idx_2 = r2_idx * merge_chunk_size + 2
            Int r2_file_idx_3 = r2_idx * merge_chunk_size + 3
            Int r2_file_idx_4 = r2_idx * merge_chunk_size + 4
            Int r2_file_idx_5 = r2_idx * merge_chunk_size + 5
            Int r2_file_idx_6 = r2_idx * merge_chunk_size + 6
            Int r2_file_idx_7 = r2_idx * merge_chunk_size + 7
            Int r2_file_idx_8 = r2_idx * merge_chunk_size + 8
            Int r2_file_idx_9 = r2_idx * merge_chunk_size + 9
            Array[File] r2_chunk = select_all([
                if (r2_file_idx_0 < round2_total_chunks) then round1_counts[r2_file_idx_0] else None,
                if (r2_file_idx_1 < round2_total_chunks) then round1_counts[r2_file_idx_1] else None,
                if (r2_file_idx_2 < round2_total_chunks) then round1_counts[r2_file_idx_2] else None,
                if (r2_file_idx_3 < round2_total_chunks) then round1_counts[r2_file_idx_3] else None,
                if (r2_file_idx_4 < round2_total_chunks) then round1_counts[r2_file_idx_4] else None,
                if (r2_file_idx_5 < round2_total_chunks) then round1_counts[r2_file_idx_5] else None,
                if (r2_file_idx_6 < round2_total_chunks) then round1_counts[r2_file_idx_6] else None,
                if (r2_file_idx_7 < round2_total_chunks) then round1_counts[r2_file_idx_7] else None,
                if (r2_file_idx_8 < round2_total_chunks) then round1_counts[r2_file_idx_8] else None,
                if (r2_file_idx_9 < round2_total_chunks) then round1_counts[r2_file_idx_9] else None
            ])

            call Merge_Whitelist_Counts as r2_merge {
                input:
                    counts_tables = r2_chunk,
                    output_path = "merged.whitelist.counts.r2." + r2_idx + ".txt.gz",
                    gzip_output = true,
                    gzip_level = 1,
                    docker = docker
            }
        }

    }
    Array[File] round2_counts = select_first([r2_merge.merged_counts, round1_counts])

    Int round3_total_chunks = length(round2_counts)
    # round 3 only runs when the previous round left more than one file
    if (round3_total_chunks > 1) {
        Int round3_num_chunks = (round3_total_chunks + merge_chunk_size - 1) / merge_chunk_size
        Array[Int] round3_indexes = range(round3_num_chunks)

        scatter (r3_idx in round3_indexes) {
            Int r3_file_idx_0 = r3_idx * merge_chunk_size + 0
            Int r3_file_idx_1 = r3_idx * merge_chunk_size + 1
            Int r3_file_idx_2 = r3_idx * merge_chunk_size + 2
            Int r3_file_idx_3 = r3_idx * merge_chunk_size + 3
            Int r3_file_idx_4 = r3_idx * merge_chunk_size + 4
            Int r3_file_idx_5 = r3_idx * merge_chunk_size + 5
            Int r3_file_idx_6 = r3_idx * merge_chunk_size + 6
            Int r3_file_idx_7 = r3_idx * merge_chunk_size + 7
            Int r3_file_idx_8 = r3_idx * merge_chunk_size + 8
            Int r3_file_idx_9 = r3_idx * merge_chunk_size + 9
            Array[File] r3_chunk = select_all([
                if (r3_file_idx_0 < round3_total_chunks) then round2_counts[r3_file_idx_0] else None,
                if (r3_file_idx_1 < round3_total_chunks) then round2_counts[r3_file_idx_1] else None,
                if (r3_file_idx_2 < round3_total_chunks) then round2_counts[r3_file_idx_2] else None,
                if (r3_file_idx_3 < round3_total_chunks) then round2_counts[r3_file_idx_3] else None,
                if (r3_file_idx_4 < round3_total_chunks) then round2_counts[r3_file_idx_4] else None,
                if (r3_file_idx_5 < round3_total_chunks) then round2_counts[r3_file_idx_5] else None,
                if (r3_file_idx_6 < round3_total_chunks) then round2_counts[r3_file_idx_6] else None,
                if (r3_file_idx_7 < round3_total_chunks) then round2_counts[r3_file_idx_7] else None,
                if (r3_file_idx_8 < round3_total_chunks) then round2_counts[r3_file_idx_8] else None,
                if (r3_file_idx_9 < round3_total_chunks) then round2_counts[r3_file_idx_9] else None
            ])

            call Merge_Whitelist_Counts as r3_merge {
                input:
                    counts_tables = r3_chunk,
                    output_path = "merged.whitelist.counts.r3." + r3_idx + ".txt.gz",
                    gzip_output = true,
                    gzip_level = 1,
                    docker = docker
            }
        }

    }
    Array[File] round3_counts = select_first([r3_merge.merged_counts, round2_counts])

    Int round4_total_chunks = length(round3_counts)
    # round 4 only runs when the previous round left more than one file
    if (round4_total_chunks > 1) {
        Int round4_num_chunks = (round4_total_chunks + merge_chunk_size - 1) / merge_chunk_size
        Array[Int] round4_indexes = range(round4_num_chunks)

        scatter (r4_idx in round4_indexes) {
            Int r4_file_idx_0 = r4_idx * merge_chunk_size + 0
            Int r4_file_idx_1 = r4_idx * merge_chunk_size + 1
            Int r4_file_idx_2 = r4_idx * merge_chunk_size + 2
            Int r4_file_idx_3 = r4_idx * merge_chunk_size + 3
            Int r4_file_idx_4 = r4_idx * merge_chunk_size + 4
            Int r4_file_idx_5 = r4_idx * merge_chunk_size + 5
            Int r4_file_idx_6 = r4_idx * merge_chunk_size + 6
            Int r4_file_idx_7 = r4_idx * merge_chunk_size + 7
            Int r4_file_idx_8 = r4_idx * merge_chunk_size + 8
            Int r4_file_idx_9 = r4_idx * merge_chunk_size + 9
            Array[File] r4_chunk = select_all([
                if (r4_file_idx_0 < round4_total_chunks) then round3_counts[r4_file_idx_0] else None,
                if (r4_file_idx_1 < round4_total_chunks) then round3_counts[r4_file_idx_1] else None,
                if (r4_file_idx_2 < round4_total_chunks) then round3_counts[r4_file_idx_2] else None,
                if (r4_file_idx_3 < round4_total_chunks) then round3_counts[r4_file_idx_3] else None,
                if (r4_file_idx_4 < round4_total_chunks) then round3_counts[r4_file_idx_4] else None,
                if (r4_file_idx_5 < round4_total_chunks) then round3_counts[r4_file_idx_5] else None,
                if (r4_file_idx_6 < round4_total_chunks) then round3_counts[r4_file_idx_6] else None,
                if (r4_file_idx_7 < round4_total_chunks) then round3_counts[r4_file_idx_7] else None,
                if (r4_file_idx_8 < round4_total_chunks) then round3_counts[r4_file_idx_8] else None,
                if (r4_file_idx_9 < round4_total_chunks) then round3_counts[r4_file_idx_9] else None
            ])

            call Merge_Whitelist_Counts as r4_merge {
                input:
                    counts_tables = r4_chunk,
                    output_path = "merged.whitelist.counts.r4." + r4_idx + ".txt.gz",
                    gzip_output = true,
                    gzip_level = 1,
                    docker = docker
            }
        }

    }
    Array[File] round4_counts = select_first([r4_merge.merged_counts, round3_counts])

    Int round5_total_chunks = length(round4_counts)
    # round 5 only runs when the previous round left more than one file
    if (round5_total_chunks > 1) {
        Int round5_num_chunks = (round5_total_chunks + merge_chunk_size - 1) / merge_chunk_size
        Array[Int] round5_indexes = range(round5_num_chunks)

        scatter (r5_idx in round5_indexes) {
            Int r5_file_idx_0 = r5_idx * merge_chunk_size + 0
            Int r5_file_idx_1 = r5_idx * merge_chunk_size + 1
            Int r5_file_idx_2 = r5_idx * merge_chunk_size + 2
            Int r5_file_idx_3 = r5_idx * merge_chunk_size + 3
            Int r5_file_idx_4 = r5_idx * merge_chunk_size + 4
            Int r5_file_idx_5 = r5_idx * merge_chunk_size + 5
            Int r5_file_idx_6 = r5_idx * merge_chunk_size + 6
            Int r5_file_idx_7 = r5_idx * merge_chunk_size + 7
            Int r5_file_idx_8 = r5_idx * merge_chunk_size + 8
            Int r5_file_idx_9 = r5_idx * merge_chunk_size + 9
            Array[File] r5_chunk = select_all([
                if (r5_file_idx_0 < round5_total_chunks) then round4_counts[r5_file_idx_0] else None,
                if (r5_file_idx_1 < round5_total_chunks) then round4_counts[r5_file_idx_1] else None,
                if (r5_file_idx_2 < round5_total_chunks) then round4_counts[r5_file_idx_2] else None,
                if (r5_file_idx_3 < round5_total_chunks) then round4_counts[r5_file_idx_3] else None,
                if (r5_file_idx_4 < round5_total_chunks) then round4_counts[r5_file_idx_4] else None,
                if (r5_file_idx_5 < round5_total_chunks) then round4_counts[r5_file_idx_5] else None,
                if (r5_file_idx_6 < round5_total_chunks) then round4_counts[r5_file_idx_6] else None,
                if (r5_file_idx_7 < round5_total_chunks) then round4_counts[r5_file_idx_7] else None,
                if (r5_file_idx_8 < round5_total_chunks) then round4_counts[r5_file_idx_8] else None,
                if (r5_file_idx_9 < round5_total_chunks) then round4_counts[r5_file_idx_9] else None
            ])

            call Merge_Whitelist_Counts as r5_merge {
                input:
                    counts_tables = r5_chunk,
                    output_path = "merged.whitelist.counts.r5." + r5_idx + ".txt.gz",
                    gzip_output = true,
                    gzip_level = 1,
                    docker = docker
            }
        }
    }
    Array[File] round5_counts = select_first([r5_merge.merged_counts, round4_counts])

    File merged_counts_output = round5_counts[0]

    scatter (i in lane_indexes) {
        call BC_Correct_Pass2 as pass2 {
            input:
                unmatched_fastq = pass1.unmatched_fastq[i],
                counts_table = merged_counts_output,
                sample_name = sample_names[i],
                whitelist = effective_whitelist,
                bc_range = effective_bc_range,
                umi_range = effective_umi_range,
                trim_len = effective_trim_len,
                min_post = min_post,
                min_trimmed_len = effective_min_trimmed_len,
                qual_tags = qual_tags,
                revcomp_trimmed = resolved_config.resolved_revcomp_trimmed,
                generate_discarded = generate_discarded,
                poly_logic = resolved_config.resolved_poly_logic,
                poly_start = resolved_config.resolved_poly_start,
                poly_window = resolved_config.resolved_poly_window,
                poly_t_min = resolved_config.resolved_poly_t_min,
                poly_a_min = resolved_config.resolved_poly_a_min,
                poly_err = resolved_config.resolved_poly_err,
                docker = docker
        }
    }

    call Summarize_Read_Accounting as read_accounting {
        input:
            sample_names = sample_names,
            trim_stats = pass1.trim_stats_line,
            pass1_read_stats = pass1.read_stats_line,
            pass2_read_stats = pass2.read_stats_line,
            docker = docker
    }

    call Plot_Read_Accounting as plot_accounting {
        input:
            per_shard = read_accounting.per_shard,
            total = read_accounting.total,
            plot_label = plot_label,
            docker = plot_docker
    }

    output {
        Array[Array[File]] fastq_pairs = transpose([pass1.matched_fastq, pass2.corrected_fastq])
        # Optional compatibility:
        Array[File] matched_fastqs = pass1.matched_fastq
        Array[File] corrected_fastqs = pass2.corrected_fastq
        Array[File] discarded_fastqs = if generate_discarded then flatten(pass2.discarded_fastqs) else []
        Array[File] unmatched_fastqs = if output_unmatched then pass1.unmatched_fastq else []
        File merged_counts_table = merged_counts_output
        File read_accounting_per_shard = read_accounting.per_shard
        File read_accounting_total = read_accounting.total
        File read_accounting_total_without_rescue = read_accounting.total_without_rescue
        File read_accounting_plots_pdf = plot_accounting.pdf
        File read_accounting_summary_png = plot_accounting.summary_png
        # The per-shard stats files the accounting is built from (same order as raw_fastqs).
        Array[File] trim_stats = pass1.trim_stats
        Array[File] pass1_read_stats = pass1.read_stats
        Array[File] pass2_read_stats = pass2.read_stats
    }
}
