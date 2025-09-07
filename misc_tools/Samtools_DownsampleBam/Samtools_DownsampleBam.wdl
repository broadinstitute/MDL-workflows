version 1.0

workflow DownsampleBam_wf {

    input {
        String sample_id
        File input_bam
        Float frac_downsample
        Int preemptible = 0

    }

    call DownsampleBam_task {
        input:
          sample_id = sample_id,
          input_bam = input_bam,
          frac_downsample = frac_downsample,
          preemptible = preemptible
    }

    output {
        File downsampled_bam = DownsampleBam_task.downsampled_bam
        File downsampled_bai = DownsampleBam_task.downsampled_bai
    }

}

task DownsampleBam_task {

    input {
        String sample_id
        File input_bam
        Float frac_downsample
        Int preemptible = 0
    }

    Int disk_space_multiplier = 2
    Int disk_space = ceil(size(input_bam, "GB")*disk_space_multiplier)
    
    command <<<
        set -ex
        
        samtools view -s ~{frac_downsample} ~{input_bam} -bo ~{sample_id}.F~{frac_downsample}.bam

        samtools index ~{sample_id}.F~{frac_downsample}.bam

    >>>

    output {
        File downsampled_bam = "~{sample_id}.F~{frac_downsample}.bam"
        File downsampled_bai = "~{sample_id}.F~{frac_downsample}.bam.bai"
    }

   runtime {
        docker:"biocontainers/samtools"
        memory: "8GB"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
    }
    
    
}
