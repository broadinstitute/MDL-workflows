version 1.0

# to make the kallisto index:

# make the cdna fasta file:
#    gffread GRCh38.gencode.v39.annotation.gtf   -g GRCh38_no_alt.fa -w GRCh38.gencode.v39.cdna.fa 
#
# build the index:
# docker run --rm \
# -u $(id -u):$(id -g) \
#  -v $(pwd):/work \
#  us-central1-docker.pkg.dev/methods-dev-lab/kallisto/kallisto \
#  kallisto index \
#    -i /work/GRCh38.gencode.v39.cdna.fa.kallisto.idx \
#    -k 31 \
#    /work/GRCh38.gencode.v39.cdna.fa



    
task KallistoQuantPairedEndSingle {
  input {
    File kallisto_index
    File r1_fastq
    File r2_fastq

    String sample_id = "sample"
    Int threads = 8
    Int bootstraps = 0               # set >0 to enable -b
    String strandedness = "unstranded"  # "unstranded" | "fr-stranded" | "rf-stranded"
    String extra_args = ""           # e.g. "--bias"
  }

  command <<<
    set -euo pipefail

    outdir="kallisto_out"
    mkdir -p "${outdir}"

    STRAND_OPT=""
    if [[ "~{strandedness}" == "fr-stranded" ]]; then
      STRAND_OPT="--fr-stranded"
    elif [[ "~{strandedness}" == "rf-stranded" ]]; then
      STRAND_OPT="--rf-stranded"
    fi

    BS_OPT=""
    if [[ "~{bootstraps}" -gt 0 ]]; then
      BS_OPT="-b ~{bootstraps}"
    fi

    kallisto quant \
      -i "~{kallisto_index}" \
      -o "${outdir}" \
      -t ~{threads} \
      ${BS_OPT} \
      ${STRAND_OPT} \
      ~{extra_args} \
      "~{r1_fastq}" "~{r2_fastq}" \
      2>&1 | tee "~{sample_id}.kallisto.log"

    # Helpful bundle for Terra download
    tar -czf "~{sample_id}.kallisto_out.tgz" -C "${outdir}" .
  >>>

  output {
    File abundance_tsv = "kallisto_out/abundance.tsv"
    File abundance_h5  = "kallisto_out/abundance.h5"
    File run_info_json = "kallisto_out/run_info.json"
    File log           = "~{sample_id}.kallisto.log"
    File output_tar    = "~{sample_id}.kallisto_out.tgz"
  }

  runtime {
    docker: "us-central1-docker.pkg.dev/methods-dev-lab/kallisto/kallisto"
    cpu: threads
    memory: "16G"
    disks: "local-disk 80 HDD"
  }
}

workflow KallistoQuantSingleSample {
  input {
    File kallisto_index
    File r1_fastq
    File r2_fastq

    String sample_id = "sample"
    Int threads = 8
    Int bootstraps = 0
    String strandedness = "unstranded"
    String extra_args = ""
  }

  call KallistoQuantPairedEndSingle {
    input:
      kallisto_index = kallisto_index,
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      sample_id = sample_id,
      threads = threads,
      bootstraps = bootstraps,
      strandedness = strandedness,
      extra_args = extra_args
  }

  output {
    File abundance_tsv = KallistoQuantPairedEndSingle.abundance_tsv
    File abundance_h5  = KallistoQuantPairedEndSingle.abundance_h5
    File run_info_json = KallistoQuantPairedEndSingle.run_info_json
    File log           = KallistoQuantPairedEndSingle.log
    File output_tar    = KallistoQuantPairedEndSingle.output_tar
  }
}
