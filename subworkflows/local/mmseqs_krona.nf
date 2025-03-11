include { MMSEQS_EASYTAXONOMY               } from '../../modules/local/mmseqs_easytaxonomy'
include { MAKE_MMSEQS_TAXONOMY_TABLE        } from '../../modules/local/make_mmseqs_taxonomy_table'
// include { KRONA_KTIMPORTTAXONOMY_MMSEQS } from '../../modules/nf-core/krona/ktimporttaxonomy/main'

workflow MMSEQS_KRONA {
    take:
    contigs   // channel: [ val(meta), path(contigs) ]
    mmseqs_db // channel: [ val(meta), path(mmseqs_db) ]
    assembler
    input_taxonomy
    // krona_db  // channel: [ val(meta), path(krona_db) ]

    // TODO: separate into subworkflows...
    main:

    ch_versions = Channel.empty()

    if (!params.taxonomy_input) {
        ch_mmseqs_files = input_taxonomy.map {}
        }
    }
    //
    // Run mmseqs easy-taxonomy on selected target database
    //
    MMSEQS_EASYTAXONOMY(
        contigs,
        mmseqs_db
    )
    ch_tophit_aln = MMSEQS_EASYTAXONOMY.out.tophit_aln.map { meta, meta2, meta3, aln ->
        def meta_new = meta + meta2 + meta3
        [meta_new, aln]
    }
    ch_tophit_report = MMSEQS_EASYTAXONOMY.out.tophit_report.map { meta, meta2, meta3, t_report ->
        def meta_new = meta + meta2 + meta3
        [meta_new, t_report]
    }
    ch_versions = ch_versions.mix(MMSEQS_EASYTAXONOMY.out.versions.first())

    ch_collected_files = ch_tophit_aln.map{ it[1] }.mix(ch_tophit_report.map{ it[1] }).collect()

    ch_mmseqs_files = ch_collected_files
        .combine(mmseqs_db)
        .map { files, mmseqs_db_map, mmseqs_db_path ->
            def meta = [mmseqs2_db_name: mmseqs_db_map.mmseqs2_db_name, assembler: assembler]
            [meta, files + [mmseqs_db_path]]
        }

    ch_mmseqs_files = ch_mmseqs_files
        .ifEmpty { 
            mmseqs_db.map { mmseqs_db_map, mmseqs_db_path ->
                [[mmseqs2_db_name: mmseqs_db_map.mmseqs2_db_name, assembler: assembler], [mmseqs_db_path]]
            }
        }

    ch_contigs_files = contigs.collect{ it[1] }.ifEmpty([])

    MAKE_MMSEQS_TAXONOMY_TABLE(
        ch_mmseqs_files,
        ch_contigs_files,
        "mmseqs_taxonomy_table.csv"
    )

    //
    // Krona plots on selected taxonomy
    //
    //TODO: add krona

    emit:
    lca           = MMSEQS_EASYTAXONOMY.out.lca           // channel: [ val(meta), [ lca ] ]
    report        = MMSEQS_EASYTAXONOMY.out.report        // channel: [ val(meta), [ report ] ]
    // krona      = KRONA_KTIMPORTTAXONOMY_
    tophit_aln    = ch_tophit_aln    // channel: [ val(meta), [ tophit_aln ] ]
    tophit_report = ch_tophit_report // channel: [ val(meta), [ tophit_report ] ]
    versions      = ch_versions                           // channel: [ versions.yml ]
}
