#!/usr/bin/env Rscript

# Load libraries
library(Seurat)
library(DropletUtils)
library(R.utils)
library(Matrix)
library(optparse)

# Parse command line arguments
option_list <- list(
  make_option(c("--indir"), type = "character", default = NULL,
              help = "Input directory containing count tarballs", metavar = "character"),
  make_option(c("--outdir"), type = "character", default = NULL,
              help = "Output directory for Seurat objects", metavar = "character"),
  make_option(c("--assign_mode"), type = "logical", default = FALSE,
              help = "Run in assign mode (requires cluster_file)", metavar = "logical"),
  make_option(c("--cluster_file"), type = "character", default = NULL,
              help = "TSV file with barcode to cluster ID assignments", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$indir)) {
  stop("--indir is required")
}
if (is.null(opt$outdir)) {
  stop("--outdir is required")
}
if (opt$assign_mode && is.null(opt$cluster_file)) {
  stop("--cluster_file is required when assign_mode is TRUE")
}

# Set paths
indir <- opt$indir
outdir <- opt$outdir
assign_mode <- opt$assign_mode
cluster_file <- opt$cluster_file

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("Starting Seurat Clustering Pipeline\n")
cat("========================================\n")
cat("Mode:", ifelse(assign_mode, "ASSIGN", "BUILD"), "\n")
cat("Input dir:", indir, "\n")
cat("Output dir:", outdir, "\n")
if (assign_mode) {
  cat("Cluster file:", cluster_file, "\n")
}
cat("========================================\n\n")

#-------------------------------
# 1. Find tarballs
#-------------------------------
tarfiles <- list.files(indir, pattern = "\\.tar\\.gz$", full.names = TRUE, recursive = TRUE)

gene_tarfiles <- tarfiles[
  grepl("gene", basename(tarfiles), ignore.case = TRUE) &
    grepl("annotated", basename(tarfiles), ignore.case = TRUE)
]
isoform_tarfiles <- tarfiles[
  grepl("isoform", basename(tarfiles), ignore.case = TRUE) &
    grepl("annotated", basename(tarfiles), ignore.case = TRUE)
]

get_sample_name <- function(path) basename(dirname(path))
gene_samples <- vapply(gene_tarfiles, get_sample_name, character(1))
isoform_samples <- vapply(isoform_tarfiles, get_sample_name, character(1))

if (length(gene_tarfiles) == 0) stop("No gene+annotated tarballs found.")
if (length(isoform_tarfiles) == 0) message("Note: No isoform+annotated tarballs found.")

cat("Found", length(gene_tarfiles), "gene tarballs\n")
cat("Found", length(isoform_tarfiles), "isoform tarballs\n\n")

#-------------------------------
# 2. Helper: reverse complement
#-------------------------------
reverse_complement <- function(seq) {
  bases <- strsplit(seq, "")[[1]]
  comp <- setNames(c("T", "G", "C", "A", "N"), c("A", "C", "G", "T", "N"))
  bases <- comp[bases]
  paste(rev(bases), collapse = "")
}

#-------------------------------
# 3. Helper: read tar, filter, return counts
#-------------------------------
read_tar_counts <- function(tarfile, sample) {
  start_time <- Sys.time()
  cat("Processing sample:", sample, "\n")
  
  outdir_tmp <- tempfile(pattern = paste0(sample, "_"))
  untar(tarfile, exdir = outdir_tmp)
  
  contents <- list.files(outdir_tmp, recursive = TRUE, full.names = TRUE)
  barcode_path <- contents[grepl("barcodes.tsv(\\.gz)?$", contents)]
  if (length(barcode_path) == 0) stop("No barcodes.tsv(.gz) found in: ", tarfile)
  data_dir <- dirname(barcode_path[1])
  
  # Ensure features.tsv.gz exists
  feat_plain <- file.path(data_dir, "features.tsv")
  feat_gz <- file.path(data_dir, "features.tsv.gz")
  genes_plain <- file.path(data_dir, "genes.tsv")
  genes_gz <- file.path(data_dir, "genes.tsv.gz")
  
  if (!file.exists(feat_gz)) {
    if (file.exists(feat_plain)) {
      R.utils::gzip(feat_plain, destname = feat_gz, overwrite = TRUE)
    } else if (file.exists(genes_gz)) {
      file.copy(genes_gz, feat_gz, overwrite = TRUE)
    } else if (file.exists(genes_plain)) {
      R.utils::gzip(genes_plain, destname = feat_gz, overwrite = TRUE)
    } else {
      stop("No features.tsv(.gz) or genes.tsv(.gz) found in: ", data_dir)
    }
  }
  
  # Read 10X
  counts <- Read10X(data.dir = data_dir)
  if (is.list(counts)) counts <- counts[[1]]
  
  # Filter cells with <100 UMIs
  totals <- Matrix::colSums(counts)
  keep <- totals >= 100
  counts <- counts[, keep, drop = FALSE]
  
  # Rename barcodes
  bc <- colnames(counts)
  bc <- sub("-[0-9]+$", "", bc)
  bc <- vapply(bc, reverse_complement, character(1))
  bc <- paste0(sample, "_", bc, "-1")
  colnames(counts) <- bc
  
  # Report
  cat("Finished sample:", sample, "| kept:", sum(keep), 
      "| filtered:", sum(!keep), "| median UMIs (kept):", 
      median(totals[keep]), "\n\n")
  
  counts
}

#-------------------------------
# 4. Build per-sample Seurat objects
#-------------------------------
build_and_save_seurat <- function(tarfiles, samples, type = "gene") {
  objs <- list()
  for (i in seq_along(tarfiles)) {
    sample <- samples[i]
    outfile <- file.path(outdir, paste0(type, "_", sample, ".rds"))
    
    if (file.exists(outfile)) {
      cat("Skipping sample", sample, "(already exists)\n")
      objs[[sample]] <- readRDS(outfile)
    } else {
      counts <- read_tar_counts(tarfiles[i], sample)
      obj <- CreateSeuratObject(counts = counts, project = sample)
      saveRDS(obj, outfile)
      objs[[sample]] <- obj
      cat("🎉 Saved Seurat object for sample:", sample,
          "| cells:", ncol(obj), "\n")
    }
  }
  return(objs)
}

gene_objs <- build_and_save_seurat(gene_tarfiles, gene_samples, "gene")
if (length(isoform_tarfiles) > 0) {
  isoform_objs <- build_and_save_seurat(isoform_tarfiles, isoform_samples, "isoform")
}

#-------------------------------
# 5. Merge objects
#-------------------------------
gene_seurat <- Reduce(function(x, y) merge(x, y), gene_objs)
if (length(isoform_tarfiles) > 0) {
  isoform_seurat <- Reduce(function(x, y) merge(x, y), isoform_objs)
}

#-------------------------------
# 6. Assign clusters if in assign mode
#-------------------------------
if (assign_mode) {
  cat("\n========================================\n")
  cat("Running in ASSIGN mode\n")
  cat("========================================\n")
  
  # Read cluster assignments
  cluster_df <- read.table(cluster_file, header = TRUE, sep = "\t", 
                           stringsAsFactors = FALSE)
  
  # Validate required columns
  if (!all(c("barcode", "cluster_id") %in% colnames(cluster_df))) {
    stop("cluster_file must contain 'barcode' and 'cluster_id' columns")
  }
  
  cat("Loaded", nrow(cluster_df), "cluster assignments\n")
  
  # Create named vector for mapping
  cluster_map <- setNames(cluster_df$cluster_id, cluster_df$barcode)
  
  # Assign clusters to gene_seurat
  barcodes <- colnames(gene_seurat)
  assigned_clusters <- cluster_map[barcodes]
  
  # Check how many matched
  n_matched <- sum(!is.na(assigned_clusters))
  cat("Matched", n_matched, "out of", length(barcodes), "barcodes\n")
  
  # Assign to metadata
  gene_seurat$assigned_cluster <- assigned_clusters
  
  # Do the same for isoform if it exists
  if (exists("isoform_seurat")) {
    barcodes_iso <- colnames(isoform_seurat)
    assigned_clusters_iso <- cluster_map[barcodes_iso]
    n_matched_iso <- sum(!is.na(assigned_clusters_iso))
    cat("Matched", n_matched_iso, "out of", length(barcodes_iso), 
        "isoform barcodes\n")
    isoform_seurat$assigned_cluster <- assigned_clusters_iso
  }
  
  # Save new assignments (barcodes that weren't in original file)
  new_barcodes <- barcodes[is.na(assigned_clusters)]
  if (length(new_barcodes) > 0) {
    new_assignments <- data.frame(
      barcode = new_barcodes,
      cluster_id = NA,
      note = "new_barcode_needs_assignment",
      stringsAsFactors = FALSE
    )
    write.table(new_assignments, 
                file.path(outdir, "new_cluster_assignments.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("Saved", length(new_barcodes), "new barcodes to new_cluster_assignments.tsv\n")
  }
}

#-------------------------------
# 7. Save merged objects
#-------------------------------
cat("\n========================================\n")
cat("Saving merged Seurat objects\n")
cat("========================================\n")

saveRDS(gene_seurat, file.path(outdir, "gene_merged_seurat.rds"))
cat("Saved gene_merged_seurat.rds | cells:", ncol(gene_seurat), "\n")

if (exists("isoform_seurat")) {
  saveRDS(isoform_seurat, file.path(outdir, "isoform_merged_seurat.rds"))
  cat("Saved isoform_merged_seurat.rds | cells:", ncol(isoform_seurat), "\n")
}

#-------------------------------
# 8. Final summary
#-------------------------------
cat("\n========================================\n")
cat("FINAL SUMMARY\n")
cat("========================================\n")
cat("Gene cells:", ncol(gene_seurat), "\n")
cat("Samples merged (gene):", length(gene_objs), "\n")
if (exists("isoform_seurat")) {
  cat("Isoform cells:", ncol(isoform_seurat), "\n")
  cat("Samples merged (isoform):", length(isoform_objs), "\n")
}
if (assign_mode) {
  cat("Mode: ASSIGN\n")
  cat("Assigned clusters added to metadata\n")
}
cat("========================================\n")
cat("Pipeline completed successfully!\n")
