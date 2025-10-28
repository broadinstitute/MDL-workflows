version 1.0

workflow SeuratClustering {
  input {
    File count_dir_tarball  # Tarball containing all count matrices
    File seurat_clustering_script  # The seurat_clustering.R script
    Boolean assign_mode = false
    File? cluster_assignment_tsv  # Required if assign_mode = true
    String docker_image = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa-seurat:latest"
    Int memory_gb = 32
    Int cpu = 4
    Int disk_gb = 100
  }

  call RunSeuratClustering {
    input:
      count_dir_tarball = count_dir_tarball,
      seurat_clustering_script = seurat_clustering_script,
      assign_mode = assign_mode,
      cluster_assignment_tsv = cluster_assignment_tsv,
      docker_image = docker_image,
      memory_gb = memory_gb,
      cpu = cpu,
      disk_gb = disk_gb
  }

  output {
    File gene_seurat_rds = RunSeuratClustering.gene_seurat_rds
    File? isoform_seurat_rds = RunSeuratClustering.isoform_seurat_rds
    File log_file = RunSeuratClustering.log_file
    File? assigned_clusters_tsv = RunSeuratClustering.assigned_clusters_tsv
  }
}

task RunSeuratClustering {
  input {
    File count_dir_tarball
    File seurat_clustering_script
    Boolean assign_mode
    File? cluster_assignment_tsv
    String docker_image
    Int memory_gb
    Int cpu
    Int disk_gb
  }

  command <<<
    set -e
    
    # Copy R script
    cp ~{seurat_clustering_script} /cromwell_root/seurat_clustering.R
    
    # Extract count directory
    mkdir -p /cromwell_root/count_data
    tar -xzf ~{count_dir_tarball} -C /cromwell_root/count_data
    
    # Copy cluster assignment file if provided
    if [ "~{assign_mode}" == "true" ]; then
      if [ -z "~{cluster_assignment_tsv}" ]; then
        echo "ERROR: assign_mode is true but no cluster_assignment_tsv provided"
        exit 1
      fi
      cp ~{cluster_assignment_tsv} /cromwell_root/cluster_assignments.tsv
    fi
    
    # Run R script
    Rscript /cromwell_root/seurat_clustering.R \
      --indir /cromwell_root/count_data \
      --outdir /cromwell_root/output \
      --assign_mode ~{assign_mode} \
      ~{if defined(cluster_assignment_tsv) then "--cluster_file /cromwell_root/cluster_assignments.tsv" else ""} \
      2>&1 | tee /cromwell_root/output/seurat_clustering.log
    
    # Ensure output files exist
    if [ ! -f /cromwell_root/output/gene_merged_seurat.rds ]; then
      echo "ERROR: gene_merged_seurat.rds was not created"
      exit 1
    fi
  >>>

  output {
    File gene_seurat_rds = "/cromwell_root/output/gene_merged_seurat.rds"
    File? isoform_seurat_rds = "/cromwell_root/output/isoform_merged_seurat.rds"
    File log_file = "/cromwell_root/output/seurat_clustering.log"
    File? assigned_clusters_tsv = if assign_mode then "/cromwell_root/output/new_cluster_assignments.tsv" else None
  }

  runtime {
    docker: docker_image
    memory: "~{memory_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_gb} SSD"
  }
}
