version 1.0

task sample_bam {
    input {
        File inputBAM
        File inputBAMindex
        String sampleName
        Float samplingRate
    }

    String output_bam_name = sample_name + ".sampled." + (if (sampling_rate > 1.0) then "estimated_" else "") + sampling_rate + ".bam"

    command <<<
        # Calculate total read count if necessary
        total_reads=$(if [[ ~{sampling_rate} > 1 ]]; then samtools view -c ~{inputBAM}; else echo 0; fi)

        # Calculate effective sampling rate
        effective_rate=$(if [[ ~{sampling_rate} > 1 ]]; then awk "BEGIN{print ~{sampling_rate} / $total_reads}"; else echo ~{sampling_rate}; fi)

        # Subsample BAM file
        samtools view -bh --subsample $effective_rate ~{inputBAM} > ~{output_bam_name}
        samtools index ~{output_bam_name}
    >>>

    output {
        File sampled_bam = output_bam_name
        File sampled_bam_index = output_bam_name + ".bai"
    }

    runtime {
        docker: "mgibio/samtools:1.16.1"
        cpu: 1
        memory: "4GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*2 + 5) + " SSD"
        preemptible: preemptible_tries
    }
}