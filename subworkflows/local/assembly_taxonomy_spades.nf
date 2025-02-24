include { SPADES as METASPADES       } from '../../modules/nf-core/spades/main'
include { SPADES as METASPADESHYBRID } from '../../modules/nf-core/spades/main'
include { MMSEQS_EASYTAXONOMY_KRONA  } from './mmseqs_easytaxonomy_krona'

workflow ASSEMBLY_TAXONOMY_SPADES {
    take:
    short_reads
    long_reads

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
            .map { id, meta_long, long_reads, meta_short, short_reads -> [meta_short, short_reads, [], long_reads] }

        METASPADESHYBRID(ch_reads_spadeshybrid, [], [])
        ch_spadeshybrid_assemblies = METASPADESHYBRID.out.scaffolds.map { meta, assembly ->
            def meta_new = meta + [assembler: "SPAdesHybrid"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spadeshybrid_assemblies)
        ch_versions = ch_versions.mix(METASPADESHYBRID.out.versions.first())
    }

    ch_mmseqs_spades_report = Channel.empty()
    if (params.spades_taxonomy) {
        MMSEQS_EASYTAXONOMY_KRONA(
            ch_assembled_contigs,
            [ [:], file(params.mmseqs2_db) ] // TODO: add meta, maybe search type or db type
        )
        ch_mmseqs_spades_report = MMSEQS_EASYTAXONOMY_KRONA.out.report
    }

    emit:
    mmseqs_report = ch_mmseqs_spades_report // channel: [ val(meta), [ mmseqs_report ] ]
    contigs       = ch_assembled_contigs    // channel: [ val(meta), [ assembly ] ]
    versions      = ch_versions             // channel: [ versions.yml ]

}
