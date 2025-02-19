include { ASSEMBLY_SPADES as MMSEQS_EASYTAXONOMY_KRONA_MEGAHIT } from '../mmseqs_easytaxonomy_krona'
include { ASSEMBLY_MEGAHIT as MMSEQS_EASYTAXONOMY_KRONA_SPADES } from '../mmseqs_easytaxonomy_krona'
include { ASSEMBLY_PREPARATION                                 } from './assembly_preparation'
include { GUNZIP as GUNZIP_ASSEMBLIES                          } from '../../modules/nf-core/gunzip'
include { QUAST                                                } from '../../modules/local/quast'

workflow CONTIGS_ASSEMBLY {
    take:
    short_read
    long_reads

    main:

    ch_versions = Channel.empty()

    // Co-assembly preparation: grouping for MEGAHIT and pooling for SPAdes
    ASSEMBLY_PREPARATION(
        short_read,
        long_reads
    )

    // Assembly

    ch_assembled_contigs = Channel.empty()

    if (!params.single_end && (!params.skip_spades || !params.skip_spadeshybrid)) {
        ASSEMBLY_SPADES(
            ASSEMBLY_PREPARATION.out.short_reads_spades,
            ASSEMBLY_PREPARATION.out.long_reads_spades
        )
        ch_assembled_contigs = ch_assembled_contigs.mix(ASSEMBLY_SPADES.out.contigs)
        ch_versions = ch_versions.mix(ASSEMBLY_SPADES.out.versions.first())
    }

    if (!params.skip_megahit) {
        ASSEMBLY_MEGAHIT(
            ASSEMBLY_PREPARATION.short_reads_grouped
        )
        ch_assembled_contigs = ch_assembled_contigs.mix(ASSEMBLY_MEGAHIT.out.contigs)
        ch_versions = ch_versions.mix(ASSEMBLY_MEGAHIT.out.versions.first())
    }

    GUNZIP_ASSEMBLIES(ch_assembled_contigs)
    ch_versions = ch_versions.mix(GUNZIP_ASSEMBLIES.out.versions)

    ch_assemblies = GUNZIP_ASSEMBLIES.out.gunzip

    ch_quast_multiqc = Channel.empty()
    if (!params.skip_quast) {
        QUAST(ch_assemblies)
        ch_versions = ch_versions.mix(QUAST.out.versions.first())

    }

    emit:
    assemblies            = ch_assemblies                      // channel: [ val(meta), [ assemblies ] ]
    quast_multiqc         = ch_quast_multiqc                   // ??? why
    spades_mmseqs_report  = ASSEMBLY_SPADES.out.mmseqs_report  // channel: [ val(meta), [ spades_report ] ]
    megahit_mmseqs_report = ASSEMBLY_MEGAHIT.out.mmseqs_report // channel: [ val(meta), [ megahit_report ] ]
    versions              = ch_versions                        // channel: [ versions.yml ]
}
