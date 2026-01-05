version 1.0

workflow STAR_align_paired_rnaseq {
  input {
    String sample_id
    
    File read1_fastq
    File read2_fastq

    # Uncompressed tar containing STAR_index_dir/ at its root
    File star_index_tar

    Int threads = 8
    String memory = "64G"
    
    # Optional: to match whatever your STAR-Fusion WDL uses
    String docker_image = "trinityctat/starfusion:latest"

    # Optional extra STAR args (ex: "--outSAMattributes NH HI AS nM --chimOutType Junctions")
    String extra_star_args = ""
  }

  call STAR_Align_SortedBam {
    input:
      sample_id = sample_id,
      read1_fastq = read1_fastq,
      read2_fastq = read2_fastq,
      star_index_tar = star_index_tar,
      threads = threads,
      memory = memory,
      docker_image = docker_image,
      extra_star_args = extra_star_args
  }

  output {
    File bam = STAR_Align_SortedBam.sorted_bam
    File bam_bai = STAR_Align_SortedBam.sorted_bam_bai
    File star_log_final = STAR_Align_SortedBam.log_final
    File star_log_out = STAR_Align_SortedBam.log_out
    File star_log_progress = STAR_Align_SortedBam.log_progress
  }
}

task STAR_Align_SortedBam {
  input {
    String sample_id
    File read1_fastq
    File read2_fastq
    File star_index_tar
    Int threads
    String memory
    String docker_image
    String extra_star_args = ""
  }

  Int disk_gb = ceil(
  (
    size(read1_fastq) +
    size(read2_fastq) +
    size(star_index_tar)
  ) / 1e9
  ) + 200
  
  command <<<
    set -euo pipefail

    mkdir -p star_index
    tar -xf "~{star_index_tar}" -C star_index

    GENOME_DIR="star_index/STAR_index_dir"
    if [[ ! -d "${GENOME_DIR}" ]]; then
      echo "ERROR: Expected STAR index directory at ${GENOME_DIR} after untarring." >&2
      echo "Contents of star_index:" >&2
      ls -lah star_index >&2 || true
      exit 1
    fi

    # Handle gzipped vs plain fastq automatically
    READ_CMD=""
    if [[ "~{read1_fastq}" == *.gz ]]; then
      READ_CMD="--readFilesCommand zcat"
    fi

    STAR \
      --runThreadN ~{threads} \
      --genomeDir "${GENOME_DIR}" \
      --readFilesIn "~{read1_fastq}" "~{read2_fastq}" \
      ${READ_CMD} \
      --outSAMtype BAM SortedByCoordinate \
      --outFileNamePrefix "~{sample_id}.star." \
      ~{extra_star_args}

    # STAR writes:
    #   star.Aligned.sortedByCoord.out.bam
    # and index may or may not be created depending on STAR version/args, so index here for consistency.
    if [[ ! -f ~{sample_id}.star.Aligned.sortedByCoord.out.bam ]]; then
      echo "ERROR: STAR did not produce ~{sample_id}.star.Aligned.sortedByCoord.out.bam" >&2
      ls -lah >&2
      exit 1
    fi

    samtools index -@ ~{threads} ~{sample_id}.star.Aligned.sortedByCoord.out.bam
  >>>

  output {
    File sorted_bam = "~{sample_id}.star.Aligned.sortedByCoord.out.bam"
    File sorted_bam_bai = "~{sample_id}.star.Aligned.sortedByCoord.out.bam.bai"
    File log_final = "~{sample_id}.star.Log.final.out"
    File log_out = "~{sample_id}.star.Log.out"
    File log_progress = "~{sample_id}.star.Log.progress.out"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: "~{memory}"
    disks: "local-disk ~{disk_gb} SSD"
  }
}
