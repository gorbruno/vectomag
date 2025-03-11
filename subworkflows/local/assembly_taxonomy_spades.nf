include { SPADES as METASPADES       } from '../../modules/nf-core/spades/main'
include { SPADES as METASPADESHYBRID } from '../../modules/nf-core/spades/main'
include { MMSEQS_KRONA  } from './mmseqs_krona'

workflow ASSEMBLY_TAXONOMY_SPADES {
    take:
    short_reads
    long_reads
    mmseqs2_db

    main:
    ch_versions = Channel.empty()
    ch_assembled_contigs = Channel.empty()

    if (!params.skip_spades) {
        METASPADES(short_reads.map { meta, reads -> [meta, reads, [], []] }, [], [])
        ch_spades_assemblies = METASPADES.out.scaffolds.map { meta, assembly ->
            def meta_new = meta + [assembler: 'SPAdes']
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spades_assemblies)
        ch_versions = ch_versions.mix(METASPADES.out.versions.first())
    }

    if (!params.skip_spadeshybrid) {
        ch_short_reads_tmp = short_reads.map { meta, reads -> [meta.id, meta, reads] }

        ch_reads_spadeshybrid = long_reads
            .map { meta, reads -> [meta.id, meta, reads] }
            .combine(ch_short_reads_tmp, by: 0)
            .map { _id, _meta_long, long_reads_in, meta_short, short_reads_in -> [meta_short, short_reads_in, [], long_reads_in] }

        METASPADESHYBRID(ch_reads_spadeshybrid, [], [])
        ch_spadeshybrid_assemblies = METASPADESHYBRID.out.scaffolds.map { meta, assembly ->
            def meta_new = meta + [assembler: "SPAdesHybrid"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spadeshybrid_assemblies)
        ch_versions = ch_versions.mix(METASPADESHYBRID.out.versions.first())
    }

    ch_mmseqs_spades_report = Channel.empty()
    if (params.spades_taxonomy && mmseqs2_db) {
        MMSEQS_KRONA(
            ch_assembled_contigs,
            mmseqs2_db,
            "SPAdes"
        )
        ch_mmseqs_spades_report = MMSEQS_KRONA.out.report
    }

    emit:
    mmseqs_report = ch_mmseqs_spades_report // channel: [ val(meta), [ mmseqs_report ] ]
    contigs       = ch_assembled_contigs    // channel: [ val(meta), [ assembly ] ]
    versions      = ch_versions             // channel: [ versions.yml ]

}
