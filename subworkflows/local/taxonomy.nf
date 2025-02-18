//
// Taxonomic sequence classification. Reads level
//

include { UNTAR as UNTAR_KRAKEN2_DB     } from '../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2               } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRONA_KRONADB                                         } from '../../modules/nf-core/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY                                } from '../../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRAKENTOOLS_KREPORT2KRONA as KREPORT2KRONA_CENTRIFUGE } from '../../modules/nf-core/krakentools/kreport2krona/main'

workflow BAM_TRIM_PRIMERS_IVAR {
    take:
    bam   // channel: [ val(meta), [ bam ], [bai] ]
    bed   // path   : bed
    fasta // channel: reference.fasta

    main:

    // Kraken2
    if (!ch_kraken2_db_file.isEmpty()) {
        if (ch_kraken2_db_file.extension in ['gz', 'tgz']) {
            // Expects to be tar.gz!
            ch_db_for_kraken2 = KRAKEN2_DB_PREPARATION(ch_kraken2_db_file).db
        }
        else if (ch_kraken2_db_file.isDirectory()) {
            ch_db_for_kraken2 = Channel
                .fromPath("${ch_kraken2_db_file}/*.k2d")
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
        }
        else {
            ch_db_for_kraken2 = Channel.empty()
        }
    }
    else {
        ch_db_for_kraken2 = Channel.empty()
    }

    KRAKEN2(
        ch_short_reads,
        ch_db_for_kraken2,
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    if ((params.centrifuge_db || params.kraken2_db) && !params.skip_krona) {
        if (params.krona_db) {
            ch_krona_db = ch_krona_db_file
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

        // Join together for Krona
        ch_tax_classifications = ch_centrifuge_for_krona
            .mix(KRAKEN2.out.results_for_krona)
            .map { classifier, meta, report ->
                def meta_new = meta + [classifier: classifier]
                [meta_new, report]
            }

        KRONA_KTIMPORTTAXONOMY(
            ch_tax_classifications,
            ch_krona_db,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
    }

    ch_versions = Channel.empty()

    //
    // iVar trim primers
    //
    IVAR_TRIM (
        bam,
        bed
    )
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())

    //
    // Count primer statistic from iVar trim
    //

    IVAR_TRIM_STATS (
        IVAR_TRIM.out.log,
        bed
    )
    ch_versions = ch_versions.mix(IVAR_TRIM_STATS.out.versions)

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        IVAR_TRIM.out.bam,
        fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = IVAR_TRIM.out.bam                    // channel: [ val(meta), bam   ]
    log_out  = IVAR_TRIM.out.log                    // channel: [ val(meta), log   ]
    primer_stats = IVAR_TRIM_STATS.out.stats        // channel: [ val(meta), stats ]
    primer_summary = IVAR_TRIM_STATS.out.summary    // channel: [ val(meta), summary ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    coverage = BAM_SORT_STATS_SAMTOOLS.out.coverage // channel: [ val(meta), [ coverage ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
