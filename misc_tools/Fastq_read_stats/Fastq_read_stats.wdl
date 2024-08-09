version 1.0

workflow Fastq_read_stats_wf {

    input {
        String sample_id
        File fastq1
        File? fastq2
        Int preemptible = 0

    }

    call Fastq_read_stats_task as fq_stats {
        input:
          sample_id = sample_id,
          fastq1 = fastq1,
          fastq2 = fastq2
    }

    output {
        File fastq_stats_file = fq_stats.fastq_stats_file
    }

}

task Fastq_read_stats_task {

    input {
        String sample_id
        File fastq1
        File? fastq2
        Int preemptible = 0
    }

    String stats_filename = "~{sample_id}.fastq.read_stats.txt"
    Int disk_space_multiplier = 3
    Int disk_space = ceil(size(fastq1, "GB")*disk_space_multiplier)
    
    command <<<

        python3 <<CODE

        import gzip, re
        sample_id = "~{sample_id}"
        fastq1 = "~{fastq1}"
        fastq2 = "~{fastq2}"
        
        filenames = [fastq1]
        if fastq2 is not None and fastq2 != "":
            filenames.append(fastq2)

        sum_seq_lens = 0
        num_seqs = 0
        num_SE_seqs = 0
        
        for fastq in filenames:
            num_SE_seqs = 0
            if re.search("\\.gz$", fastq):
                fh = gzip.open(fastq, 'rt')
            else:
                fh = open(fastq, 'rt')
            linecounter = 0
            for line in fh:
                linecounter += 1
                if linecounter % 4 == 2:
                    # seqline
                    line = line.rstrip()
                    sum_seq_lens += len(line)
                    num_seqs += 1
                    num_SE_seqs += 1
        
        mean_read_length = sum_seq_lens / num_seqs


        
        with open("~{stats_filename}", "wt") as ofh:
            print("num_SE_seqs:\t{}\tsum_bases:\t{}\tmean_seq_len:\t{}".format(num_SE_seqs,
                                                                               sum_seq_lens,
                                                                               round(sum_seq_lens/num_seqs)), file=ofh)

        CODE

    >>>

    output {
        File fastq_stats_file = "~{stats_filename}"
    }

   runtime {
        docker:"ubuntu:20.04"
        memory: "8GB"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
    
    
}

