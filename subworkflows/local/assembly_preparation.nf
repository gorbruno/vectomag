include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS } from '../../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                            } from '../../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS         } from '../../modules/local/pool_single_reads'

workflow ASSEMBLY_PREPARATION {
    take:
    short_reads     // channel: [ val(meta), path(reads1), path(reads2) ]
    long_reads

    main:

    ch_versions = Channel.empty()

    // Co-assembly preparation: grouping for MEGAHIT and pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_grouped = short_reads
            .map { meta, reads -> [meta.group, meta, reads] }
            .groupTuple(by: 0)
            .map { group, _metas, reads ->
                def assemble_as_single = params.single_end || (params.bbnorm && params.coassemble_group)
                def meta = [:]
                meta.id = "group-${group}"
                meta.group = group
                meta.single_end = assemble_as_single
                if (assemble_as_single) {
                    [meta, reads.collect { it }, []]
                }
                else {
                    [meta, reads.collect { it[0] }, reads.collect { it[1] }]
                }
            }
        // long reads
        // group and set group as new id
        ch_long_reads_grouped = long_reads
            .map { meta, reads -> [meta.group, meta, reads] }
            .groupTuple(by: 0)
            .map { group, _metas, reads ->
                def meta = [:]
                meta.id = "group-${group}"
                meta.group = group
                [meta, reads.collect { it }]
            }
    }
    else {
        ch_short_reads_grouped = short_reads
            .filter { it[0].single_end }
            .map { meta, reads -> [meta, [reads], []] }
            .mix(
                short_reads.filter { !it[0].single_end }.map { meta, reads -> [meta, [reads[0]], [reads[1]]] }
            )
        ch_long_reads_grouped = long_reads
    }

    if (!params.skip_spades || !params.skip_spadeshybrid) {
        if (params.coassemble_group) {
            if (params.bbnorm) {
                ch_short_reads_spades = ch_short_reads_grouped.map { [it[0], it[1]] }
            }
            else {
                POOL_SHORT_SINGLE_READS(
                    ch_short_reads_grouped.filter { it[0].single_end }
                )
                POOL_PAIRED_READS(
                    ch_short_reads_grouped.filter { !it[0].single_end }
                )
                ch_short_reads_spades = POOL_SHORT_SINGLE_READS.out.reads.mix(POOL_PAIRED_READS.out.reads)
            }
        }
        else {
            ch_short_reads_spades = short_reads
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid) {
            POOL_LONG_READS(ch_long_reads_grouped)
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        }
        else {
            ch_long_reads_spades = Channel.empty()
        }
    }
    else {
        ch_short_reads_spades = Channel.empty()
        ch_long_reads_spades = Channel.empty()
    }

    emit:
    short_reads_spades  = ch_short_reads_spades  // channel: [ val(meta), [ lca ] ]
    long_reads_spades   = ch_long_reads_spades   // channel: [ val(meta), [ report ] ]
    short_reads_grouped = ch_short_reads_grouped // channel: [ val(meta), [ tophit_aln ] ]
    versions            = ch_versions            // channel: [ versions.yml ]
}
