include { MEGAHIT                   } from '../../modules/nf-core/megahit/main'
include { MMSEQS_EASYTAXONOMY_KRONA } from './mmseqs_easytaxonomy_krona'

workflow ASSEMBLY_MEGAHIT {
    take:
    short_reads

    main:
    ch_versions = Channel.empty()
    ch_assembled_contigs = Channel.empty()

    MEGAHIT(short_reads)

    ch_megahit_assemblies = MEGAHIT.out.contigs.map { meta, assembly ->
        def meta_new = meta + [assembler: 'MEGAHIT']
        [meta_new, assembly]
    }
    ch_assembled_contigs = ch_assembled_contigs.mix(ch_megahit_assemblies)

    ch_mmseqs_megahit_report = Channel.empty()
    if (params.megahit_taxonomy) {
      MMSEQS_EASYTAXONOMY_KRONA(
        ch_assembled_contigs,
        [ [:], file(params.mmseqs2_db) ] // TODO: add meta, maybe search type or db type
      )
      ch_mmseqs_megahit_report = MMSEQS_EASYTAXONOMY_KRONA.out.report
    }

    emit:
    mmseqs_report = ch_mmseqs_spades     // channel: [ val(meta), [ mmseqs_report ] ]
    contigs       = ch_assembled_contigs // channel: [ val(meta), [ assembly ] ]
    versions      = ch_versions          // channel: [ versions.yml ]

}
