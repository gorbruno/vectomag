include { ASSEMBLY_SPADES               } from '../assembly_spades'
include { ASSEMBLY_MEGAHIT              } from '../assembly_megahit'
// include { KRONA_KTIMPORTTAXONOMY_MMSEQS } from '../../modules/nf-core/krona/ktimporttaxonomy/main'

workflow ASSEMBLY {
    take:
    contigs   // channel: [ val(meta), path(contigs) ]
    mmseqs_db // channel: [ val(meta), path(mmseqs_db) ]
    // krona_db  // channel: [ val(meta), path(krona_db) ]

    // TODO: separate into subworkflows...
    main:

    ch_versions = Channel.empty()

    //
    // Run mmseqs easy-taxonomy on selected target database
    //
    MMSEQS_EASYTAXONOMY(
        contigs,
        mmseqs_db
    )
    ch_versions = ch_versions.mix(MMSEQS_EASYTAXONOMY.out.versions.first())

    //
    // Krona plots on selected taxonomy
    //
    //TODO: test and add krona

    emit:
    lca           = MMSEQS_EASYTAXONOMY.out.lca           // channel: [ val(meta), [ lca ] ]
    report        = MMSEQS_EASYTAXONOMY.out.report        // channel: [ val(meta), [ report ] ]
    // krona      = KRONA_KTIMPORTTAXONOMY_
    tophit_aln    = MMSEQS_EASYTAXONOMY.out.tophit_aln    // channel: [ val(meta), [ tophit_aln ] ]
    tophit_report = MMSEQS_EASYTAXONOMY.out.tophit_report // channel: [ val(meta), [ tophit_report ] ]
    versions      = ch_versions                           // channel: [ versions.yml ]
}
