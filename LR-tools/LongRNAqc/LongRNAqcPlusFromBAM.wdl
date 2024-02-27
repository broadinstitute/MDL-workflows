version 1.0


import "../utilities/sampleBAM.wdl" as sampleBAM
import "../sqanti3/sqanti3FromBAM.wdl" as sqanti3FromBAMWorkflow
import "LongRNAqcFromBAM.wdl" as LongRNAqcFromBAMWorkflow
import "../IsoQuant/IsoQuantMakeDB.wdl" as IsoQuantMakeDBWorkflow
import "../IsoQuant/IsoQuantQuantify.wdl" as IsoQuantQuantifyWorkflow


struct SampleBamAndIndex {
    String sample_name
    File bam
    File bam_index
}



workflow LongRNAqcPlusFromBam {

    meta {
        description: "Using a BAM as input, run the pacbio adapted version of rnaseqc, Sqanti3, and IsoQuant."
    }

    input {
        Array[String] sampleName
        String dataType
        String ?strandedness
        Array[File] inputBAM
        Array[File] inputBAMIndex
        String chromosomesList # comma seprarated
        File referenceFasta
        File referenceGTF
        File ?referenceGTF_DB
        File collapsedReferenceGTF
        File cagePeak
        File polyAMotifs
        Float ?samplingRate
        String BAMToGTFConversionMethod
        String transcriptQuantification = "with_ambiguous"
        String geneQuantification = "with_inconsistent"
        Boolean allowNonPrimary
        Int cpu = 4
        Int preemptible_tries = 3

    }

    scatter(i in range(length(sampleName))) {
        SampleBamAndIndex sampleBamAndIndex = object { 
            sample_name: sampleName[i], 
            bam: inputBAM[i],
            bam_index: inputBAMIndex[i]
        }
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

    # scatter(sample in createStructTask.sampleBamAndIndex) {
    scatter(sample in sampleBamAndIndex) {
        if (defined(samplingRate)) {
            Float sampling_rate = select_first([samplingRate])
            call sampleBAM.sample_bam as sampleBam {
                input:
                    sampleName = sample.sample_name,
                    inputBAM = sample.bam,
                    inputBAMindex = sample.bam_index,
                    samplingRate = sampling_rate,
                    maxRetries = preemptible_tries
            }
        }

        File bam_file = select_first([sampleBam.sampled_bam, sample.bam])
        File bam_file_index = select_first([sampleBam.sampled_bam_index, sample.bam_index])

        call LongRNAqcFromBAMWorkflow.LongRNAqc as LongRNAqc {
            input:
                sampleName = sample.sample_name,
                inputBAM = bam_file,
                inputBAMIndex = bam_file_index,
                collapsedReferenceGTF = collapsedReferenceGTF,
                maxRetries = preemptible_tries
        }

        call sqanti3FromBAMWorkflow.sqanti3FromBam as sqanti3FromBam {
            input:
                sampleName = sample.sample_name,
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

        call IsoQuantQuantifyWorkflow.isoquantQuantify as isoquantQuantify {
            input:
                sampleName = sample.sample_name,
                inputBAM = bam_file,
                inputBAMIndex = bam_file_index,
                referenceFasta = referenceFasta,
                referenceAnnotation = isoquantDB,
                dataType = dataType,
                strandedness = strandedness,
                transcriptQuantification = transcriptQuantification,
                geneQuantification = geneQuantification,
                noModelConstruction = true,
                preemptible_tries = preemptible_tries
        }
    }

    # Here is where the plotting would happen

    output {
        # Array[SampleBamAndIndex] sampleBamAndIndex = createStructTask.sampleBamAndIndex
        Array[File] sampledBAM = bam_file
        Array[File] sampledBAMindex= bam_file_index
        Array[File] rnaseqc_gene_reads_gct = LongRNAqc.rnaseqc_gene_reads_gct
        Array[File] rnaseqc_gene_fragments_gct = LongRNAqc.rnaseqc_gene_fragments_gct
        Array[File] rnaseqc_gene_tpm_gct = LongRNAqc.rnaseqc_gene_tpm_gct
        Array[File] rnaseqc_exon_reads_gct = LongRNAqc.rnaseqc_exon_reads_gct
        Array[File] rnaseqc_exon_cv_tsv = LongRNAqc.rnaseqc_exon_cv_tsv
        Array[File] rnaseqc_metrics_tsv = LongRNAqc.rnaseqc_metrics_tsv
        Array[File] sqantiClassificationTSV = sqanti3FromBam.sqantiClassificationTSV
        Array[File] sqantiJunctionsTSV = sqanti3FromBam.sqantiJunctionsTSV
        Array[File] sqantiCorrectedFasta = sqanti3FromBam.sqantiCorrectedFasta
        Array[File] sqantiCorrectedGTF = sqanti3FromBam.sqantiCorrectedGTF
        Array[Array[File]] isoquantOutputs = isoquantQuantify.isoquantOutputs
        Array[File?] isoquantTranscriptModelsGTF = isoquantQuantify.transcriptModelsGTF
        Array[File?] isoquantReadAssignmentsTSV = isoquantQuantify.readAssignmentsTSV
    }

}




