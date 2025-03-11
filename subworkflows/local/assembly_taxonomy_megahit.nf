include { MEGAHIT                   } from '../../modules/nf-core/megahit/main'
include { MMSEQS_KRONA } from './mmseqs_krona'

workflow ASSEMBLY_TAXONOMY_MEGAHIT {
    take:
    short_reads
    mmseqs2_db

    main:
    ch_versions = Channel.empty()
    ch_assembled_contigs = Channel.empty()

    MEGAHIT(short_reads)

    ch_contigs = MEGAHIT.out.contigs


    ch_megahit_assemblies = ch_contigs.map { meta, assembly ->
        def meta_new = meta + [assembler: 'MEGAHIT']
        [meta_new, assembly]
    }

    ch_assembled_contigs = ch_assembled_contigs.mix(ch_megahit_assemblies)

    ch_mmseqs_megahit_report = Channel.empty()
    if (params.megahit_taxonomy && mmseqs2_db) {
        MMSEQS_KRONA(
            ch_assembled_contigs,
            mmseqs2_db,
            "MEGAHIT"
        )
        ch_mmseqs_megahit_report = MMSEQS_KRONA.out.report
    }


    emit:
    mmseqs_report = ch_mmseqs_megahit_report // channel: [ val(meta), [ mmseqs_report ] ]
    contigs       = ch_assembled_contigs     // channel: [ val(meta), [ assembly ] ]
    versions      = ch_versions              // channel: [ versions.yml ]

}
