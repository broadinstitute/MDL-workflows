version 1.0

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
        Boolean qual_tags = false
        Boolean revcomp_trimmed = false
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
    String trim_poly_logic_arg = if (poly_enabled) then ("-l " + resolved_poly_logic) else ""
    String poly_start_arg = if (poly_start != "") then ("--poly-start " + poly_start) else ""
    String poly_window_arg = if (poly_window != "") then ("--poly-window " + poly_window) else ""
    String poly_t_min_arg = if (poly_t_min != "") then ("--polyT-min " + poly_t_min) else ""
    String poly_a_min_arg = if (poly_a_min != "") then ("--polyA-min " + poly_a_min) else ""
    String poly_err_arg = if (poly_err != "") then ("--poly-err " + poly_err) else ""
    String poly_logic_arg = if (poly_enabled) then ("--poly-logic " + resolved_poly_logic) else ""
    String qual_tags_arg = if (qual_tags) then "--qual-tags" else ""
    String revcomp_arg = if (revcomp_trimmed) then "--revcomp-trimmed" else ""

    String trimmed_prefix = sample_name
    String counts_path = trimmed_prefix + ".whitelist.counts.txt"
    String matched_fastq = trimmed_prefix + ".matched.fastq.gz"
    String unmatched_fastq = trimmed_prefix + ".unmatched.fastq.gz"
    String counts_table_name = counts_path

    Int diskGB = select_first([disk_size_gb, ceil(size(input_fastq, "GB") * 3 + 20)])

    command <<<
        set -euo pipefail

        if [[ -z "~{whitelist}" ]]; then
          echo "Error: whitelist is required. Provide library_config (resolved whitelist) or whitelist input." >&2
          exit 1
        fi

        /opt/pipeline/fastq_preprocessing_pipeline/trim_split_nocap_extracttags.sh \
          --stdout-found \
          -i ~{input_fastq} \
          -o ~{trimmed_prefix} \
          -t 1 \
          -a ~{end_adapter} \
          -g ~{start_adapter} \
          -e ~{endedness} \
          ~{trim_poly_logic_arg} \
        | bc_correct_split \
          --fastq4 \
          --whitelist ~{whitelist} \
          --outdir . \
          --threads 1 \
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
    >>>

    output {
        File matched_fastq = "~{matched_fastq}"
        File unmatched_fastq = "~{unmatched_fastq}"
        File counts_table = "~{counts_table_name}"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "4 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Merge_Whitelist_Counts {
    input {
        Array[File] counts_tables
        String output_path = "merged.whitelist.counts.txt"
        String docker
        String? countsum_script
        Int? disk_size_gb
    }

    String effective_script = select_first([countsum_script, "/opt/pipeline/fastq_preprocessing_pipeline/sum_whitelist_counts.py"])
    Int diskGB = select_first([disk_size_gb, ceil(size(counts_tables, "GB") * 3 + 20)])

    command <<<
        set -euo pipefail
        python3 ~{effective_script} --output ~{output_path} ~{sep=" " counts_tables}
    >>>

    output {
        File merged_counts = "~{output_path}"
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


task Resolve_BC_Params_From_Config {
    input {
        String? library_config
        File? library_config_file
        String docker
    }

    command <<<
        set -euo pipefail

        has_preset="~{if (defined(library_config)) then "1" else "0"}"
        has_file="~{if (defined(library_config_file)) then "1" else "0"}"
        if [ "$has_preset" = "1" ] && [ "$has_file" = "1" ]; then
          echo "Error: provide exactly one of library_config (preset) or library_config_file." >&2
          exit 1
        fi
        if [ "$has_preset" = "0" ] && [ "$has_file" = "0" ]; then
          echo "Error: library_config (preset) or library_config_file is required." >&2
          exit 1
        fi

        if [ "$has_file" = "1" ]; then
          config_path="~{library_config_file}"
        else
          preset_name="~{library_config}"
          preset_base="${preset_name##*/}"
          preset_base="${preset_base%.yaml}"
          preset_base="${preset_base%.yml}"
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

        python3 - <<'PY'
        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "pipeline_wrapper", "/opt/pipeline/fastq_preprocessing_pipeline/pipeline_wrapper.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        config = module.load_config(config_path)
        params = module.build_from_config(config)

        def write(path, value, default="50"):
            with open(path, "w", encoding="utf-8") as f:
                if value is None:
                    f.write(str(default))
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
        write("resolved_poly_logic.txt", params.get("poly", {}).get("logic"), default="")
        write("resolved_poly_start.txt", params.get("poly", {}).get("start"), default="")
        write("resolved_poly_window.txt", params.get("poly", {}).get("window"), default="")
        write("resolved_poly_t_min.txt", params.get("poly", {}).get("min_t"), default="")
        write("resolved_poly_a_min.txt", params.get("poly", {}).get("min_a"), default="")
        write("resolved_poly_err.txt", params.get("poly", {}).get("err"), default="")
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
        Boolean qual_tags = false
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
    String discard_arg = if (generate_discarded) then "--discarded" else ""

    String output_prefix = sample_name
    String corrected_fastq = output_prefix + ".corrected.fastq.gz"

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
          --threads 1 \
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
    >>>

    output {
        File corrected_fastq = "~{corrected_fastq}"
        Array[File] discarded_fastqs = glob("~{output_prefix}*.discarded*.fastq.gz")
    }

    runtime {
        docker: docker
        cpu: 1
        memory: "4 GB"
        disks: "local-disk ~{diskGB} SSD"
        preemptible: 3
    }
}


workflow BC_Barcode_Extract_And_Correct_Array {
    meta {
        description: "Extract and correct cell barcodes/UMIs in two passes over arrays of FASTQs."
    }

    input {
        Array[File] raw_fastqs
        Array[String] sample_names
        File? whitelist
        Float min_post = 0.975
        String? library_config
        File? library_config_file

        Boolean qual_tags = false
        Boolean revcomp_trimmed = false
        Boolean generate_discarded = false
        Boolean output_unmatched = false

        File? whitelist_override

        String docker

    }

    Array[Int] lane_indexes = range(length(raw_fastqs))
    call Resolve_BC_Params_From_Config as resolved_config {
        input:
            library_config = library_config,
            library_config_file = library_config_file,
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
                input_fastq = raw_fastqs[i],
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
                revcomp_trimmed = revcomp_trimmed,
                poly_logic = resolved_config.resolved_poly_logic,
                poly_start = resolved_config.resolved_poly_start,
                poly_window = resolved_config.resolved_poly_window,
                poly_t_min = resolved_config.resolved_poly_t_min,
                poly_a_min = resolved_config.resolved_poly_a_min,
                poly_err = resolved_config.resolved_poly_err,
                docker = docker
        }
    }

    call Merge_Whitelist_Counts as merged_counts {
        input:
            counts_tables = pass1.counts_table,
            docker = docker,
            output_path = "merged.whitelist.counts.txt"
    }

    scatter (i in lane_indexes) {
        call BC_Correct_Pass2 as pass2 {
            input:
                unmatched_fastq = pass1.unmatched_fastq[i],
                counts_table = merged_counts.merged_counts,
                sample_name = sample_names[i],
                whitelist = effective_whitelist,
                bc_range = effective_bc_range,
                umi_range = effective_umi_range,
                trim_len = effective_trim_len,
                min_post = min_post,
                min_trimmed_len = effective_min_trimmed_len,
                qual_tags = qual_tags,
                revcomp_trimmed = revcomp_trimmed,
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

    output {
        Array[File] matched_fastqs = pass1.matched_fastq
        Array[File] corrected_fastqs = pass2.corrected_fastq
        Array[File] discarded_fastqs = if generate_discarded then flatten(pass2.discarded_fastqs) else []
        Array[File] unmatched_fastqs = if output_unmatched then pass1.unmatched_fastq else []
        File merged_counts_table = merged_counts.merged_counts
    }
}
