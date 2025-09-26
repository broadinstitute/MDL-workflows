WDL added here partition_bam_by_id.wdl which is a task within the main LRAA wdl that calls the script:
https://github.com/MethodsDev/LongReadAlignmentAssembler/blob/main/util/sc/partition_bam_by_cell_cluster.py
Originally created by bhaas@broadinstitute.org, the script above partitions a single cell Kinnex pre-processed bam, into individual clusters when provided with barcode to id mapping.

We're re-purposing it here to partition the bam by patient id for BrCA.
