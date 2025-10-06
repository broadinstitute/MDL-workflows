# Evaluating transcript alignment coverage by reads across the transcript lengths.

- run minimap2 to align the reads to your cDNA fasta sequences

 minimap2 -a -N0 -t 10 ref_annot.cdna.fa reads.fastq.gz | samtools view -bo cdna_aligned.mm2.bam


- get transcript alignment depth profiles.

  samtools depth cdna_aligned.mm2.bam > cdna_aligned.mm2.bam.depth

- extract coverage profiles

   eval_isoform_alignment_coverage_profile.py -b cdna_aligned.mm2.bam \
                                              -d cdna_aligned.mm2.bam.depth \
                                              -o cdna_aligned.mm2.bam.depth.cov_info

    
    