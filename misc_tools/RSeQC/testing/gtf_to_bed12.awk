BEGIN {
  FS = "\t";
  OFS = "\t";
}

$0 !~ /^#/ && $3 == "exon" {
  transcript_id = "";
  gene_id = "";

  if (match($9, /transcript_id "([^"]+)"/, transcript_match)) {
    transcript_id = transcript_match[1];
  }
  if (match($9, /gene_id "([^"]+)"/, gene_match)) {
    gene_id = gene_match[1];
  }
  if (transcript_id == "") {
    next;
  }

  key = transcript_id;
  chrom[key] = $1;
  strand[key] = $7;
  name[key] = gene_id "|" transcript_id;

  start = $4 - 1;
  end = $5;

  if (!(key in tx_start) || start < tx_start[key]) {
    tx_start[key] = start;
  }
  if (!(key in tx_end) || end > tx_end[key]) {
    tx_end[key] = end;
  }

  exon_count[key]++;
  exon_starts[key, exon_count[key]] = start;
  exon_ends[key, exon_count[key]] = end;
  transcripts[key] = 1;
}

END {
  for (key in transcripts) {
    n = exon_count[key];

    for (i = 1; i <= n; i++) {
      for (j = i + 1; j <= n; j++) {
        if (exon_starts[key, j] < exon_starts[key, i]) {
          tmp_start = exon_starts[key, i];
          tmp_end = exon_ends[key, i];
          exon_starts[key, i] = exon_starts[key, j];
          exon_ends[key, i] = exon_ends[key, j];
          exon_starts[key, j] = tmp_start;
          exon_ends[key, j] = tmp_end;
        }
      }
    }

    block_sizes = "";
    block_starts = "";
    for (i = 1; i <= n; i++) {
      block_sizes = block_sizes (exon_ends[key, i] - exon_starts[key, i]) ",";
      block_starts = block_starts (exon_starts[key, i] - tx_start[key]) ",";
    }

    print chrom[key], tx_start[key], tx_end[key], name[key], 0, strand[key], \
      tx_start[key], tx_end[key], 0, n, block_sizes, block_starts;
  }
}
