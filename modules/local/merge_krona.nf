process MERGE_KRONA {
    label 'process_single'

    conda "bioconda::krona=2.8.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1':
        'biocontainers/krona:2.8.1--pl5321hdfd78af_1' }"

    input:
    val report_prefix
    path taxonomy, stageAs: 'taxonomy.tab'

    output:
    path "sum.read.taxonomy.krona.html", emit: html
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '2.8.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    TAXONOMY=\$(find -L . -name '*.tab' -exec dirname {} \\;)
    echo \$TAXONOMY

    ktImportTaxonomy \\
        $args \\
        -o sum.read.taxonomy.krona.html \\
        -tax \$TAXONOMY/ \\
        $report_prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: $VERSION
    END_VERSIONS
    """
}
