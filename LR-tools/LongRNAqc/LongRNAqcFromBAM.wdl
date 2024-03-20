version 1.0


task LongRNAqcTask {
    input {
        String sampleName
        File collapsedReferenceGTF
        File inputBAM
        File inputBAMIndex

        String docker
        Int cpu
        String memoryGB
        Int preemptible
        Int maxRetries
        Float diskSpaceMultiplier
        File monitoringScript = "gs://mdl-refs/util/cromwell_monitoring_script2.sh"
    }


    Int disk_space = ceil( (size(inputBAM, "GB") + size(collapsedReferenceGTF, "GB") ) * diskSpaceMultiplier)

    command <<<
        bash ~{monitoringScript} > monitoring.log &

        # Check if reference_annotation is provided as a .gz file
        ref_annotation_arg="~{collapsedReferenceGTF}"
        if [[ "~{collapsedReferenceGTF}" == *.gz ]]; then
            # Extract the .gz file and update reference_annotation variable
            gunzip -c "~{collapsedReferenceGTF}" > temp_reference_annotation.gtf
            ref_annotation_arg="temp_reference_annotation.gtf"
        fi


        ln -s ~{inputBAM} ~{sampleName}.bam
        ln -s ~{inputBAMIndex} ~{sampleName}.bam.bai

        rnaseqc ${ref_annotation_arg} ~{sampleName}.bam . -u


        mv ~{sampleName}.bam.gene_reads.gct ~{sampleName}.rnaseqc.gene_reads.gct
        mv ~{sampleName}.bam.gene_fragments.gct ~{sampleName}.rnaseqc.gene_fragments.gct
        mv ~{sampleName}.bam.gene_tpm.gct ~{sampleName}.rnaseqc.gene_tpm.gct
        mv ~{sampleName}.bam.exon_reads.gct ~{sampleName}.rnaseqc.exon_reads.gct
        mv ~{sampleName}.bam.exon_cv.tsv ~{sampleName}.rnaseqc.exon_cv.tsv
        mv ~{sampleName}.bam.metrics.tsv ~{sampleName}.rnaseqc.metrics.tsv

        gzip ~{sampleName}.rnaseqc.*
    >>>

    output {
        File rnaseqc_gene_reads_gct = "~{sampleName}.rnaseqc.gene_reads.gct.gz"
        File rnaseqc_gene_fragments_gct = "~{sampleName}.rnaseqc.gene_fragments.gct.gz"
        File rnaseqc_gene_tpm_gct = "~{sampleName}.rnaseqc.gene_tpm.gct.gz"
        File rnaseqc_exon_reads_gct = "~{sampleName}.rnaseqc.exon_reads.gct.gz"
        File rnaseqc_exon_cv_tsv = "~{sampleName}.rnaseqc.exon_cv.tsv.gz"
        File rnaseqc_metrics_tsv = "~{sampleName}.rnaseqc.metrics.tsv.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker}"
        disks: "local-disk " + disk_space + " HDD"
        memory: "~{memoryGB} GiB"
        cpu: cpu
        preemptible: preemptible
        maxRetries: maxRetries
    }
}
   
workflow LongRNAqc {

    input {
        String sampleName
        File inputBAM
        File inputBAMIndex
        File collapsedReferenceGTF

        String docker = "trinityctat/rnaseqc_plus:latest"
        Int cpu = 10
        String memoryGB="50"
        Int preemptible = 0
        Int maxRetries = 0
        Float diskSpaceMultiplier = 4.0
    }

    call LongRNAqcTask {
        input:
           sampleName = sampleName,
           collapsedReferenceGTF = collapsedReferenceGTF,
           inputBAM = inputBAM,
           inputBAMIndex = inputBAMIndex,

           docker = docker,
           cpu = 1,
           memoryGB = "10",
           preemptible = preemptible,
           maxRetries = maxRetries,
           diskSpaceMultiplier = diskSpaceMultiplier 

    }

    output {
        File rnaseqc_gene_reads_gct = LongRNAqcTask.rnaseqc_gene_reads_gct
        File rnaseqc_gene_fragments_gct = LongRNAqcTask.rnaseqc_gene_fragments_gct
        File rnaseqc_gene_tpm_gct = LongRNAqcTask.rnaseqc_gene_tpm_gct
        File rnaseqc_exon_reads_gct = LongRNAqcTask.rnaseqc_exon_reads_gct
        File rnaseqc_exon_cv_tsv = LongRNAqcTask.rnaseqc_exon_cv_tsv
        File rnaseqc_metrics_tsv = LongRNAqcTask.rnaseqc_metrics_tsv
    }
}
 