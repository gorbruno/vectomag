process MMSEQS_EASYTAXONOMY {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::mmseqs2=17.b804f"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1'
        : 'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(db_target)

    output:
    tuple val(meta), path("${prefix}_lca.tsv")      , emit: lca
    tuple val(meta), path("${prefix}_report")       , emit: report
    tuple val(meta), path("${prefix}_tophit_aln")   , emit: tophit_aln
    tuple val(meta), path("${prefix}_tophit_report"), emit: tophit_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    # Extract files with specified args based suffix | remove suffix | isolate longest common substring of files
    DB_TARGET_PATH_NAME=\$(find -L "${db_target}/" -maxdepth 1 -name "${args2}" | sed 's/\\.[^.]*\$//' | sed -e 'N;s/^\\(.*\\).*\\n\\1.*\$/\\1\\n\\1/;D' )

    mmseqs \\
        easy-taxonomy \\
        ${fasta} \\
        \$DB_TARGET_PATH_NAME \\
        ${prefix} \\
        tmp1 \\
        ${args} \\
        --threads ${task.cpus}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: "*.dbtype"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: /')
    END_VERSIONS
    """
}
