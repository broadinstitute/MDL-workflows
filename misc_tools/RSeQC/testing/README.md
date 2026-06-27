# RSeQC gene body coverage test

This test uses the small SIRV BAM and annotation from:

```text
~/projects/LRAA_PAPER_Analyses/iso-reconstruct-benchmark2/workflows/StringTie_docker/testing/data/
```

The test data files are copied into `RSeQC/testing/data/`. The RSeQC gene-body coverage module requires a 12-column BED gene model, so `gtf_to_bed12.awk` converts the exon-only SIRV GTF into transcript-level BED12.

To regenerate the BED12 file:

```bash
awk -f RSeQC/testing/gtf_to_bed12.awk \
  RSeQC/testing/data/SIRVs1-7.annot.reduced.gtf \
  | sort -k1,1 -k2,2n > RSeQC/testing/data/SIRVs1-7.annot.reduced.bed12
```

To run the test:

```bash
RSeQC/testing/run_test.sh
```
