version 1.0


import "../sqanti3/sqanti3FromBAM.wdl" as sqanti3FromBAMWorkflow
import "LongRNAqcFromBAM.wdl" as LongRNAqcFromBAMWorkflow
import "../IsoQuant/IsoQuantMakeDB.wdl" as IsoQuantMakeDBWorkflow
import "../IsoQuant/IsoQuantQuantify.wdl" as IsoQuantQuantifyWorkflow


struct SampleBamAndIndex {
    String sample_name
    File bam
    File bam_index
}



workflow rnaseqcPlusFromBam {

    meta {
        description: "Using a BAM as input, run the pacbio adapted version of rnaseqc, Sqanti3, and IsoQuant."
    }

    input {
        Array[String] sampleName
        String dataType
        Array[File] inputBAM
        Array[File] inputBAMIndex
        String chromosomesList # comma seprarated
        File referenceFasta
        File referenceGTF
        File ?referenceGTF_DB
        File collapsedReferenceGTF
        File cagePeak
        File polyAMotifs
        String BAMToGTFConversionMethod
        Boolean allowNonPrimary
        Int cpu = 4

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
                    isCompleteGeneDB = false
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
                isCompleteGeneDB = true
        }
    }

    File isoquantDB = select_first([isoquantMakeGeneDB_fromRef.geneDB, isoquantMakeGeneDB_fromDB.geneDB, referenceGTF_DB])

    # scatter(sample in createStructTask.sampleBamAndIndex) {
    scatter(sample in sampleBamAndIndex) {
        call LongRNAqcFromBAMWorkflow.LongRNAqc as LongRNAqc {
            input:
                sampleName = sample.sample_name,
                inputBAM = sample.bam,
                inputBAMIndex = sample.bam_index,
                collapsedReferenceGTF = collapsedReferenceGTF
        }

        call sqanti3FromBAMWorkflow.sqanti3FromBam as sqanti3FromBam {
            input:
                sampleName = sample.sample_name,
                inputBAM = sample.bam,
                inputBAMIndex = sample.bam_index,
                chromosomesList = chromosomesList,
                referenceGTF = referenceGTF,
                referenceFasta = referenceFasta,
                cagePeak = cagePeak,
                polyAMotifs = polyAMotifs,
                conversionMethod = BAMToGTFConversionMethod,
                allowNonPrimary = allowNonPrimary
        }

        call IsoQuantQuantifyWorkflow.isoquantQuantify as isoquantQuantify {
            input:
                sampleName = sample.sample_name,
                inputBAM = sample.bam,
                inputBAMIndex = sample.bam_index,
                referenceFasta = referenceFasta,
                geneDB = isoquantDB,
                dataType = dataType,
                noModelConstruction = true
        }
    }

    # Here is where the plotting would happen

    output {
        # Array[SampleBamAndIndex] sampleBamAndIndex = createStructTask.sampleBamAndIndex
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
        Array[File] isoquantReadAssignmentsTSV = isoquantQuantify.readAssignmentsTSV
    }

}




