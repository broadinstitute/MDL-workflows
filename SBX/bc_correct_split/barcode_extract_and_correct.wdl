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
    String counts_path = trimmed_prefix + ".whitelist.counts.txt.gz"
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
        predefinedMachineType: "n2d-custom-2-4096"
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
        cpu: 1
        memory: "4 GB"
        predefinedMachineType: "n2d-custom-1-4096"
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
        predefinedMachineType: "n2d-custom-1-4096"
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
                revcomp_trimmed = resolved_config.resolved_revcomp_trimmed,
                poly_logic = resolved_config.resolved_poly_logic,
                poly_start = resolved_config.resolved_poly_start,
                poly_window = resolved_config.resolved_poly_window,
                poly_t_min = resolved_config.resolved_poly_t_min,
                poly_a_min = resolved_config.resolved_poly_a_min,
                poly_err = resolved_config.resolved_poly_err,
                docker = docker
        }
    }

    Int merge_chunk_size = 10
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

    Array[File] round2_counts = r2_merge.merged_counts

    Int round3_total_chunks = length(round2_counts)
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

    Array[File] round3_counts = r3_merge.merged_counts

    Int round4_total_chunks = length(round3_counts)
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

    Array[File] round4_counts = r4_merge.merged_counts

    Int round5_total_chunks = length(round4_counts)
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

    File merged_counts_output = r5_merge.merged_counts[0]

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

    output {
        Array[Array[File]] fastq_pairs = transpose([pass1.matched_fastq, pass2.corrected_fastq])
        # Optional compatibility:
        Array[File] matched_fastqs = pass1.matched_fastq
        Array[File] corrected_fastqs = pass2.corrected_fastq
        Array[File] discarded_fastqs = if generate_discarded then flatten(pass2.discarded_fastqs) else []
        Array[File] unmatched_fastqs = if output_unmatched then pass1.unmatched_fastq else []
        File merged_counts_table = merged_counts_output
    }
}
