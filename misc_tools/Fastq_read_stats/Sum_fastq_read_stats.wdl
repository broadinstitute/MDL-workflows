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

    command <<<

        python3 - ~{sep=" " fastq_stats_files} <<CODE
        import sys

        total_reads = 0

        for stats_file in sys.argv[1:]:
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
