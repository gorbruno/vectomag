include { ASSEMBLY_TAXONOMY_SPADES    } from './assembly_taxonomy_spades'
include { ASSEMBLY_TAXONOMY_MEGAHIT   } from './assembly_taxonomy_megahit'
include { ASSEMBLY_PREPARATION        } from './assembly_preparation'
include { GUNZIP as GUNZIP_ASSEMBLIES } from '../../modules/nf-core/gunzip'
include { QUAST                       } from '../../modules/local/quast'

workflow CONTIGS_ASSEMBLY {
    take:
    short_reads
    long_reads
    mmseqs2_db

    main:

    ch_versions = Channel.empty()

    // Co-assembly preparation: grouping for MEGAHIT and pooling for SPAdes
    ASSEMBLY_PREPARATION(
        short_reads,
        long_reads
    )

    // Assembly

    ch_assembled_contigs = Channel.empty()

    ch_spades_mmseqs_report = Channel.empty()
    if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)) {
        ASSEMBLY_TAXONOMY_SPADES(
            ASSEMBLY_PREPARATION.out.short_reads_spades,
            ASSEMBLY_PREPARATION.out.long_reads_spades,
            mmseqs2_db
        )
        ch_assembled_contigs = ch_assembled_contigs.mix(ASSEMBLY_TAXONOMY_SPADES.out.contigs)
        ch_spades_mmseqs_report = ASSEMBLY_TAXONOMY_SPADES.out.mmseqs_report
        ch_versions = ch_versions.mix(ASSEMBLY_TAXONOMY_SPADES.out.versions.first())
    }

    ch_megahit_mmseqs_report = Channel.empty()
    if (!params.skip_megahit) {
        ASSEMBLY_TAXONOMY_MEGAHIT(
            ASSEMBLY_PREPARATION.out.short_reads_grouped,
            mmseqs2_db
        )
        ch_assembled_contigs = ch_assembled_contigs.mix(ASSEMBLY_TAXONOMY_MEGAHIT.out.contigs)
        ch_megahit_mmseqs_report = ASSEMBLY_TAXONOMY_MEGAHIT.out.mmseqs_report
        ch_versions = ch_versions.mix(ASSEMBLY_TAXONOMY_MEGAHIT.out.versions.first())
    }

    GUNZIP_ASSEMBLIES(ch_assembled_contigs)
    ch_versions = ch_versions.mix(GUNZIP_ASSEMBLIES.out.versions)

    ch_assemblies = GUNZIP_ASSEMBLIES.out.gunzip

    ch_quast_report = Channel.empty()
    if (!params.skip_quast) {
        QUAST(ch_assemblies)
        ch_quast_report = QUAST.out.report
        ch_versions = ch_versions.mix(QUAST.out.versions.first())

    }

    emit:
    assemblies            = ch_assemblies            // channel: [ val(meta), [ assemblies ] ]
    quast_report          = ch_quast_report          // channel: [ val(meta), [ quast_report ] ]
    spades_mmseqs_report  = ch_spades_mmseqs_report  // channel: [ val(meta), [ spades_report ] ]
    megahit_mmseqs_report = ch_megahit_mmseqs_report // channel: [ val(meta), [ megahit_report ] ]
    versions              = ch_versions              // channel: [ versions.yml ]
}
