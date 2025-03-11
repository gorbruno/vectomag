/*
 * geNomad: Identification of mobile genetic elements
 */

include { GENOMAD_DOWNLOAD   } from '../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND   } from '../../modules/nf-core/genomad/endtoend/main'

workflow VIRUS_IDENTIFICATION {
    take:
    ch_assemblies   // [ [ meta] , fasta    ], input scaffolds (mandatory)
    ch_genomad_db   // [ db                 ], presupplied geNomad database

    main:
    ch_versions = Channel.empty()

    ch_identified_viruses = Channel.empty()
    if (!ch_genomad_db.isEmpty()) {
        ch_identified_viruses = GENOMAD_ENDTOEND ( ch_assemblies, ch_genomad_db ).virus_fasta
        ch_versions.mix( GENOMAD_ENDTOEND.out.versions )
    }

    emit:
    identified_viruses = ch_identified_viruses
    versions = ch_versions

}
