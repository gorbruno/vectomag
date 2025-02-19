//
// Taxonomic sequence classification. Reads level
//

include { UNTAR as UNTAR_KRAKEN2_DB                             } from '../../modules/nf-core/untar/main'
include { UNTAR as CENTRIFUGEDB_UNTAR                           } from '../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2                            } from '../../modules/nf-core/kraken2/kraken2/main'
include { CENTRIFUGE_CENTRIFUGE as CENTRIFUGE                   } from '../../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT                                    } from '../../modules/nf-core/centrifuge/kreport/main'
include { KRONA_KRONADB                                         } from '../../modules/nf-core/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY                                } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { MERGE_KRONA                                           } from '../../modules/local/merge_krona'
include { KRAKENTOOLS_KREPORT2KRONA as KREPORT2KRONA_KRAKEN     } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENTOOLS_KREPORT2KRONA as KREPORT2KRONA_CENTRIFUGE } from '../../modules/nf-core/krakentools/kreport2krona/main'

// TODO: separate subworkflows for centrifuge and kraken2

workflow READ_TAXONOMY {
    take:
    kraken2_db  // path: kraken2_db
    krona_db    // path: krona_db
    short_reads // channel: [ val(meta), path(reads) ]

    main:

    ch_versions = Channel.empty()

    // Centrifuge
    if (!params.centrifuge_db) {
        ch_db_for_centrifuge = Channel.empty()
    }
    else {
        if (file(params.centrifuge_db).isDirectory()) {
            ch_db_for_centrifuge = Channel.of(file(params.centrifuge_db, checkIfExists: true))
        }
        else {
            ch_db_for_centrifuge = CENTRIFUGEDB_UNTAR(Channel.of([[id: 'db'], file(params.centrifuge_db, checkIfExists: true)])).untar.map { it[1] }.first()
            ch_versions = ch_versions.mix(CENTRIFUGEDB_UNTAR.out.versions.first())
        }
    }

    CENTRIFUGE(
        short_reads,
        ch_db_for_centrifuge,
        false,
        false,
    )
    ch_versions = ch_versions.mix(CENTRIFUGE.out.versions.first())

    CENTRIFUGE_KREPORT(CENTRIFUGE.out.results, ch_db_for_centrifuge)
    ch_versions = ch_versions.mix(CENTRIFUGE_KREPORT.out.versions.first())

    // Kraken2
    if (!kraken2_db.isEmpty()) {
        if (kraken2_db.extension in ['gz', 'tgz']) {
            // Expects to be tar.gz!
            UNTAR_KRAKEN2_DB([ kraken2_db.simpleName, kraken2_db ])
            ch_db_for_kraken2 = UNTAR_KRAKEN2_DB.out.untar
            ch_versions = ch_versions.mix(UNTAR_KRAKEN2_DB.out.versions)
        }
        else if (kraken2_db.isDirectory()) {

            /*
            TODO: rewrite?
            ch_db_for_kraken2 = Channel
                .fromPath("${kraken2_db}/*.k2d")
                .collect()
                .map { file ->
                    if (file.size() >= 3) {
                        def db_name = file[0].getParent().getName()
                        [db_name, file]
                    }
                    else {
                        error("Kraken2 requires '{hash,opts,taxo}.k2d' files.")
                    }
                }
            */

            Channel
                .fromPath("${kraken2_db}/*.k2d")
                .collect()
                .map { file ->
                    if (file.size() < 3) {
                        error("Kraken2 requires '{hash,opts,taxo}.k2d' files.")
                    }
                }
            ch_db_for_kraken2 = kraken2_db

        }
        else {
            ch_db_for_kraken2 = Channel.empty()
        }
    }
    else {
        ch_db_for_kraken2 = Channel.empty()
    }

    KRAKEN2(
        short_reads,
        ch_db_for_kraken2,
        false,
        false,
    )
    ch_kraken2_for_krona = KRAKEN2.out.report.map { meta, files -> ['kraken2', meta, files] }
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())



    if ((params.centrifuge_db || params.kraken2_db) && !params.skip_krona) {
        if (params.krona_db) {
            ch_krona_db = krona_db
        }
        else {
            KRONA_KRONADB()
            ch_krona_db = KRONA_KRONADB.out.db
            ch_versions = ch_versions.mix(KRONA_KRONADB.out.versions)
        }

        if (params.centrifuge_db) {
            ch_centrifuge_for_krona = KREPORT2KRONA_CENTRIFUGE(CENTRIFUGE_KREPORT.out.kreport).txt.map { meta, files -> ['centrifuge', meta, files] }
            ch_versions = ch_versions.mix(KREPORT2KRONA_CENTRIFUGE.out.versions.first())
        }
        else {
            ch_centrifuge_for_krona = Channel.empty()
        }

        /*
        TODO: add option to krona from txt?
        if (params.kraken_txt) {
            ch_kraken2_for_krona = ch_kraken2_for_krona
                .mix(KREPORT2KRONA_KRAKEN(KRAKEN2.out.report).txt.map { meta, files -> ['kraken2', meta, files] })
            ch_versions = ch_versions.mix(KREPORT2KRONA_KRAKEN.out.versions.first())
        }
        */

        // Join together for Krona
        ch_tax_classifications = ch_centrifuge_for_krona
            .mix(ch_kraken2_for_krona)
            .map { classifier, meta, report ->
                def meta_new = meta + [classifier: classifier]
                [meta_new, report]
            }

        if (params.krona_sum) {
            ch_krona_reports = ch_tax_classifications
                .map { meta, report ->
                    def sample_name = meta.id
                    "${report},${sample_name}"
                }
                .collect()
                .map { it.join(' ') } // TODO: add meta info about classifier

            MERGE_KRONA(
                ch_krona_reports,
                ch_krona_db
            )
            ch_versions = ch_versions.mix(MERGE_KRONA.out.versions)
        }
        else {
            KRONA_KTIMPORTTAXONOMY(
                ch_tax_classifications,
                ch_krona_db,
            )
            ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
        }
    }

    emit:
    kraken2_report = KRAKEN2.out.report                 // channel: [ val(meta) report   ]
    centrifuge_report  = CENTRIFUGE_KREPORT.out.kreport // channel: [ val(meta), log   ]

    versions = ch_versions                              // channel: [ versions.yml ]
}
