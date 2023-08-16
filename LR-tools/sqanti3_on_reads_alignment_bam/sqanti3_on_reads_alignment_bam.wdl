version 1.0


task split_bam_per_chromosome {
    input {
        File inputBAM
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    command <<<
        mkdir -p split_dir
        split_bam_per_chromosome.sh ~{inputBAM} split_dir
    >>>

    output {
        Array[File] chromosomeBAMs = glob("split_dir/*.bam")
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*2 + 10) + " HDD"
        docker: docker
    }
}


task backformatBAM {
    input {
        File inputBAM
        String outputType = "sam"       # sam or bam
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    command <<<
        baseBamName=$(basename ~{inputBAM} | sed 's/\(.*\)\..*/\1/')
        echo "${baseBamName}.~{outputType}" > output_name.txt

        reformat.sh \
            in=~{inputBAM} \
            out=${baseBamName}.~{outputType} \
            sam=1.3
    >>>

    output {
        File backformatedBAM = read_string("output_name.txt")
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputBAM, "GB")*10 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_CTATLR {
    input {
        File inputSAM
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    String alignmentGTF_name = basename("~{inputSAM}", ".sam")

    command <<<
        baseSamName=$(basename ~{inputSAM} | sed 's/\(.*\)\..*/\1/')

        SAM_to_gxf.pl --format gtf --allow_non_primary \
            --sam ~{inputSAM} \
            > temp.sam_to_gxf
        grep -P "[A-z0-1]" temp.sam_to_gxf > ${baseSamName}.gtf
    >>>

    output {
        File alignmentGTF = "~{alignmentGTF_name}"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputSAM, "GB")*3 + 10) + " HDD"
        docker: docker
    }
}


task convertSAMtoGTF_cDNACupcake {
    input {
        File inputSAM
        File? reference_fasta
        Boolean correct_fasta = false
        Int memoryGB
        # Int diskSizeGB
        String docker
    }

    String extra_arg = if correct_fasta then "--reference_genome ~{reference_fasta} --fasta_correction" else ""
    String alignmentGTF_name = basename("~{inputSAM}", ".sam")

    command <<<
        baseSamName=$(basename ~{inputSAM} | sed 's/\(.*\)\..*/\1/')

        convert_SAM_to_GTF_for_SQANTI3.py \
            --sam_file ~{inputSAM} \
            --output_prefix ${baseSamName} \
            ~{extra_arg}
    >>>

    output {
        File alignmentGTF = "~{alignmentGTF_name}"
        File? correctedFasta = glob("*corrected.fasta")
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk " + ceil(size(inputSAM, "GB")*3 + 10) + " HDD"
        docker: docker
    }
}


task concatenate_gtfs {
    input {
        Array[File] files
        Int memoryGB
        # Int diskSizeGB
        # String docker
    }

    command <<<
        cat '~{sep="' '" files}' > "concatenated.gtf"
    >>>

    output {
        File concatenatedGTF = "concatenated.gtf"
    }

    runtime {
        cpu: 1
        memory: "~{memoryGB} GiB"
        disks: "local-disk 100 HDD"
        docker: "alpine:latest"
    }
}


task run_sqanti {
    input {
        File input_gtf
        File reference_gtf
        File reference_fasta
        File cage_peak
        File polyA_motifs
        Int cpu
        Int memoryGB
        Int diskSizeGB
        String docker
    }

    command <<<
        sqanti3_qc.py \
            --report both \
            --chunks ${cpu}
            --dir sqanti_out_dir \
            --CAGE_peak ${cage_peak} \
            --polyA_motif_list ${polyA_motifs} \
            --skipORF \
            --window 20 \
            --isoform_hits \
            ${input_gtf} \
            ${reference_gtf} \
            ${reference_fasta}
    >>>

    output {
        Array[File] sqanti_outputs = glob("sqanti_out_dir/*")
    }

    runtime {
        cpu: cpu
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
        docker: docker
    }
}


workflow sqanti3_on_reads_alignment_bam {

    meta {
        description: "Run Sqanti3 classification of non-assembled reads on an input BAM with the alignments without having to rerun the alignment."
    }

    input {
        File inputBAM
        String conversion_method = "CTAT-LR"
        File reference_gtf
        File reference_fasta
        File cage_peak
        File polyA_motifs
        Int cpu = 4
        Int memoryGB = 16
        Int diskSizeGB = 500
    }

    String docker = "us-east4-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus@sha256:7f04c4ffef87d29f198442ab8da80cd1eea9679b54884a4bf21a6a123c0bf856"
    File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

    call split_bam_per_chromosome {
        input:
            inputBAM = inputBAM,
            memoryGB = memoryGB,
            # diskSizeGB = ,
            docker = docker
    }

    scatter(chromosomeBAM in split_bam_per_chromosome.chromosomeBAMs) {
        call backformatBAM {
            input:
                inputBAM = chromosomeBAM,
                memoryGB = memoryGB,
                # diskSizeGB = diskSizeGB,
                docker = docker
        }

        if (conversion_method == "CTAT-LR") {
            call convertSAMtoGTF_CTATLR {
                input:
                    inputSAM = backformatBAM.backformatedBAM,
                    memoryGB = memoryGB,
                    # diskSizeGB = diskSizeGB,
                    docker = docker
            }
        }

        if (conversion_method == "cDNACupcake") {
            call convertSAMtoGTF_cDNACupcake {
                input:
                    inputSAM = backformatBAM.backformatedBAM,
                    reference_fasta = reference_fasta,
                    memoryGB = memoryGB,
                    # diskSizeGB = diskSizeGB,
                    docker = docker
            }
        }

        File converted_gtf = select_first([convertSAMtoGTF_CTATLR.alignmentGTF, convertSAMtoGTF_cDNACupcake.alignmentGTF])
    }

    call concatenate_gtfs {
        input:
            files = converted_gtf,
            memoryGB = 8
    }

    call run_sqanti {
        input:
            input_gtf = concatenate_gtfs.concatenatedGTF,
            reference_gtf = reference_gtf,
            reference_fasta = reference_fasta,
            cage_peak = cage_peak,
            polyA_motifs = polyA_motifs,
            cpu = cpu,
            memoryGB = 32,
            diskSizeGB = diskSizeGB,
            docker = docker
    }

    output {
        Array[File] sqanti_outputs = run_sqanti.sqanti_outputs
    }
}

