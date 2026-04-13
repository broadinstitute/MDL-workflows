version 1.0

workflow Sum_fastq_read_stats_wf {

    input {
        Array[File] fastq_stats_files
        Int preemptible = 0
    }

    call Sum_fastq_read_stats_task {
        input:
          fastq_stats_files = fastq_stats_files,
          preemptible = preemptible
    }

    output {
        Int total_reads = Sum_fastq_read_stats_task.total_reads
    }

}

task Sum_fastq_read_stats_task {

    input {
        Array[File] fastq_stats_files
        Int preemptible = 0
    }

    File fastq_stats_list = write_lines(fastq_stats_files)

    command <<<

        python3 <<CODE
        total_reads = 0

        with open("~{fastq_stats_list}", "rt") as list_fh:
            stats_files = [line.strip() for line in list_fh if line.strip()]

        for stats_file in stats_files:
            with open(stats_file, "rt") as fh:
                fields = fh.read().strip().split("\t")

            try:
                field_index = fields.index("num_SE_seqs:")
            except ValueError:
                raise RuntimeError(f"num_SE_seqs field not found in {stats_file}")

            if field_index + 1 >= len(fields):
                raise RuntimeError(f"num_SE_seqs value missing in {stats_file}")

            total_reads += int(fields[field_index + 1])

        with open("total_reads.txt", "wt") as ofh:
            print(total_reads, file=ofh)
        CODE

    >>>

    output {
        Int total_reads = read_int("total_reads.txt")
    }

   runtime {
        docker:"python"
        memory: "2 GB"
        predefinedMachineType: "n2d-highcpu-2"
        bootDiskSizeGb: 12
        disks: "local-disk 10 HDD"
        cpu: 2
        preemptible: preemptible
    }

}
