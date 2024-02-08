version 1.0


task stringtie_merge {
    input {
        Array[File] gtf_assemblies
        File ?reference_annotation
        Int cpu
        Int memoryGB
        Int numThreads
        Int diskSizeGB
        String docker
    }

    command <<<
        # Initialize variable for reference annotation argument
        ref_annotation_arg=""

        # Check if reference_annotation is provided and is a .gz file
        if [ -n "~{reference_annotation}" ]; then
            if [[ "~{reference_annotation}" == *.gz ]]; then
                # Extract the .gz file and update reference_annotation variable
                gunzip -c "~{reference_annotation}" > temp_reference_annotation.gtf
                ref_annotation_arg="-G temp_reference_annotation.gtf"
            else
                # Use the provided reference_annotation file directly
                ${ref_annotation_arg}="-G ~{reference_annotation}"
            fi
        fi

        # Generate a new list of filenames, decompressing if necessary
        for gtf in '~{sep="' '" gtf_assemblies}'; do
            if [[ "$gtf" == *.gz ]]; then
                gunzip -c "$gtf" > "${gtf%.gz}"
                echo "${gtf%.gz}" >> assembly_GTF_list.txt
            else
                echo "$gtf" >> assembly_GTF_list.txt
            fi
        done

        stringtie \
            --merge \
            -p ~{numThreads} \
            -i \
            -o stringtie_merged.gtf \
            ${ref_annotation_arg} $(cat assembly_GTF_list.txt)
    >>>

    output {
        File merged_assembly = "stringtie_merged.gtf"
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }

}

task stringtie_re_estimate {
    input {
        String sample_name
        File bam_alignment      # single file because called through scatter()
        String reads_type
        File merged_assembly
        Boolean outputCoverage
        Boolean outputGeneAbundances
        Int cpu
        Int memoryGB
        Int numThreads
        Int diskSizeGB
        String docker
    }

    String outputName = "StringTie_out_reestimated_~{sample_name}.gtf"
    String ballgown_out = if outputCoverage then "-B" else ""
    String gene_abundance = if outputGeneAbundances then "-A StringTie_out_reestimated_~{sample_name}.gene_abund.tab" else ""

    command <<<
        stringtie \
            -e ~{ballgown_out} \
            -G ~{merged_assembly} \
            -p ~{numThreads} \
            -o ~{outputName} ~{gene_abundance} \
            ~{reads_type} \
            ~{bam_alignment}
    >>>

    output {
        File reestimated_assembly = "~{outputName}"
        File? coverageOutput = "~{sample_name}.gene_abund.tab"
        Array[File]? ctabOutputs =  glob("*.ctab")
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow stringtie_merge {

    meta {
        description: "Merge all the given StringTie assemblies into a unified assembly, then re-estimate transcript abundance for each run using the bam alignement that produced the assembly."
    }

    input {
        File ?reference_annotation       # this.samples.reference_annotation or any of them, which should all be the same, ideally as a workspace.referenceData_mm38_ field
        Array[String] sample_names      # this.samples.sample_id
        Array[File] gtf_assemblies      # this.samples.stringTieGTF
        Array[File] bam_alignments      # this.samples.input_bam
        Boolean outputCoverage
        Boolean outputGeneAbundances
        # ideally these 2 Boolean could be Arrays[Boolean] in case of different types of runs being merged, but maybe a single Array[String] with the option to use directly would work too, using Array requires nested zip() calls for scatter() however
        Boolean has_longreads               # true/false
        Boolean has_shortreads              # true/false
        Int cpus = 4
        Int memoryGB = 16
        Int diskSizeGB = 500
    }

    String docker = "us.gcr.io/broad-dsde-methods/kockan/stringtie@sha256:ca2a163c7acdcacba741ea98d573080c15f153de18bd1566d24e8d2f1729ce89"
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    # combine all StringTie assemblies into a unique one
    call stringtie_merge {
        input:
            gtf_assemblies = gtf_assemblies,
            reference_annotation = reference_annotation,
            numThreads = cpus * 2,
            cpu = cpus,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB,
            docker = docker
    }

    #if has_longreads && has_shortreads then "--mix"
    #if has_longreads && !has_shortreads then "-L"
    # if !has_longreads && has_shortreads then default
    String reads_type = if has_longreads then (if has_shortreads then "--mix" else "-L") else ""

    # may need to have a map sample_name to bam_alignment? to scatter on a single object
    scatter(pair in zip(sample_names, bam_alignments)) {
        call stringtie_re_estimate {
            input:
                sample_name = pair.left,
                bam_alignment = pair.right,
                reads_type = reads_type,
                merged_assembly = stringtie_merge.merged_assembly,
                outputCoverage = outputCoverage,
                outputGeneAbundances = outputGeneAbundances,
                numThreads = cpus*2,
                cpu = cpus,
                memoryGB = memoryGB,
                diskSizeGB = diskSizeGB,
                docker = docker
        }
    }

    output {
        File merged_assembly = stringtie_merge.merged_assembly
        Array[File] reestimated_assemblies = stringtie_re_estimate.reestimated_assembly
        Array[File?] coverageOutput = stringtie_re_estimate.coverageOutput
        Array[Array[File]?] ctabOutputs = stringtie_re_estimate.ctabOutputs
    }
}

