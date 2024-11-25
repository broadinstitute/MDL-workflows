version 1.0


import "../utilities/sampleBAM.wdl" as sampleBAM
import "../sqanti3/sqanti3FromBAM.wdl" as sqanti3FromBAMWorkflow
import "LongRNAqcFromBAM.wdl" as LongRNAqcFromBAMWorkflow
import "../IsoQuant/IsoQuantMakeDB.wdl" as IsoQuantMakeDBWorkflow
import "../IsoQuant/IsoQuantQuantify.wdl" as IsoQuantQuantifyWorkflow


workflow LongRNAqcPlusFromBam {

    meta {
        description: "Using a BAM as input, run the pacbio adapted version of rnaseqc, Sqanti3, and IsoQuant."
    }

    input {
        String sampleName
        String dataType
        String ?strandedness
        File inputBAM
        File inputBAMIndex
        String chromosomesList # comma seprarated
        File referenceFasta
        File referenceGTF
        File ?referenceGTF_DB
        File collapsedReferenceGTF
        File ?cagePeak
        File ?polyAMotifs
        Float ?samplingRate
        String BAMToGTFConversionMethod
        String transcriptQuantification = "unique_only"
        String geneQuantification = "unique_splicing_consistent"
        Boolean allowNonPrimary
        Boolean runLongRNAqc = true
        Boolean runSqanti = true
        Boolean runIsoQuant = false
        Int preemptible_tries = 3
    }


    if (defined(referenceGTF_DB)) {
        String db_filename = basename(select_first([referenceGTF_DB]))

        # If file has extension 'gtf', remove the extension such that the returned string does not match `db_filename`
        if (sub(db_filename, "gtf", "") != db_filename) {
            call IsoQuantMakeDBWorkflow.isoquantMakeGeneDB as isoquantMakeGeneDB_fromDB {
                input:
                    gtfToDB = select_first([referenceGTF_DB]),
                    isCompleteGeneDB = false,
                    preemptible_tries = preemptible_tries
            }
        }
        #if (sub(db_filename, "gtf", "") == db_filename) {
        #    File providedIsoquantDB = referenceGTF_DB
        #}
    }
    if (!defined(referenceGTF_DB)) {
        call IsoQuantMakeDBWorkflow.isoquantMakeGeneDB as isoquantMakeGeneDB_fromRef {
            input:
                gtfToDB = referenceGTF,
                isCompleteGeneDB = true,
                preemptible_tries = preemptible_tries
        }
    }

    File isoquantDB = select_first([isoquantMakeGeneDB_fromRef.geneDB, isoquantMakeGeneDB_fromDB.geneDB, referenceGTF_DB])

### might need to prefilter based on allowNonPrimary in case they can get sampled


    if (defined(samplingRate)) {
        Float sampling_rate = select_first([samplingRate])
        call sampleBAM.sample_bam as sampleBam {
            input:
                sampleName = sampleName,
                inputBAM = inputBAM,
                inputBAMindex = inputBAMIndex,
                samplingRate = sampling_rate,
                maxRetries = preemptible_tries
        }
    }

    File bam_file = select_first([sampleBam.sampled_bam, inputBAM])
    File bam_file_index = select_first([sampleBam.sampled_bam_index, inputBAM])

    if (runLongRNAqc) {
        call LongRNAqcFromBAMWorkflow.LongRNAqc as LongRNAqc {
            input:
                sampleName = sampleName,
                inputBAM = bam_file,
                inputBAMIndex = bam_file_index,
                collapsedReferenceGTF = collapsedReferenceGTF,
                maxRetries = preemptible_tries
        }
    }

    if (runSqanti) {
        call sqanti3FromBAMWorkflow.sqanti3FromBam as sqanti3FromBam {
            input:
                sampleName = sampleName,
                inputBAM = bam_file,
                inputBAMIndex = bam_file_index,
                chromosomesList = chromosomesList,
                referenceGTF = referenceGTF,
                referenceFasta = referenceFasta,
                cagePeak = cagePeak,
                polyAMotifs = polyAMotifs,
                conversionMethod = BAMToGTFConversionMethod,
                allowNonPrimary = allowNonPrimary,
                preemptible_tries = preemptible_tries
        }
    }

    if (runIsoQuant) {
        call IsoQuantQuantifyWorkflow.isoquantQuantify as isoquantQuantify {
            input:
                sampleName = sampleName,
                inputBAM = bam_file,
                inputBAMIndex = bam_file_index,
                referenceFasta = referenceFasta,
                referenceAnnotation = isoquantDB,
                dataType = dataType,
                strandedness = strandedness,
                transcriptQuantification = transcriptQuantification,
                geneQuantification = geneQuantification,
                noModelConstruction = false,
                preemptible_tries = preemptible_tries
        }
    }


    # Here is where the plotting would happen

    output {
        File sampledBAM = bam_file
        File sampledBAMindex = bam_file_index
        File ?rnaseqc_gene_reads_gct = LongRNAqc.rnaseqc_gene_reads_gct
        File ?rnaseqc_gene_fragments_gct = LongRNAqc.rnaseqc_gene_fragments_gct
        File ?rnaseqc_gene_tpm_gct = LongRNAqc.rnaseqc_gene_tpm_gct
        File ?rnaseqc_exon_reads_gct = LongRNAqc.rnaseqc_exon_reads_gct
        File ?rnaseqc_exon_cv_tsv = LongRNAqc.rnaseqc_exon_cv_tsv
        File ?rnaseqc_metrics_tsv = LongRNAqc.rnaseqc_metrics_tsv
        File ?sqantiClassificationTSV = sqanti3FromBam.sqantiClassificationTSV
        File ?sqantiJunctionsTSV = sqanti3FromBam.sqantiJunctionsTSV
        File ?sqantiCorrectedFasta = sqanti3FromBam.sqantiCorrectedFasta
        File ?sqantiCorrectedGTF = sqanti3FromBam.sqantiCorrectedGTF
        Array[File] ?allIsoquantOutputs = isoquantQuantify.allIsoquantOutputs
        File ?referenceCountsTSV = isoquantQuantify.referenceCountsTSV
        File ?referenceReadAssignmentsTSV = isoquantQuantify.referenceReadAssignmentsTSV
        File ?constructedTranscriptModelsGTF = isoquantQuantify.constructedTranscriptModelsGTF
        File ?constructedTranscriptCountsTSV = isoquantQuantify.constructedTranscriptCountsTSV
        File ?constructedTranscriptReadAssignmentsTSV = isoquantQuantify.constructedTranscriptReadAssignmentsTSV
        File ?groupedReferenceGeneCountsTSV = isoquantQuantify.groupedReferenceGeneCountsTSV
        File ?groupedReferenceTranscriptCountsTSV = isoquantQuantify.groupedReferenceTranscriptCountsTSV
        File ?groupedConstructedTranscriptCountsTSV = isoquantQuantify.groupedConstructedTranscriptCountsTSV
    }

}




