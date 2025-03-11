process QUAST {
    tag "${meta.assembler}-${meta.id}"
    label "process_medium"

    conda "bioconda::quast=5.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(assembly)

    output:
    path "QUAST/*"                       , emit: qc
    path "QUAST/report_rawassemblies.tsv", emit: report
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args   ?: ''
    """
    metaquast.py \\
        --threads "${task.cpus}" \\
        --rna-finding \\
        --max-ref-number 0 \\
        -l "${meta.assembler}-${meta.id}" \\
        "${assembly}" \\
        -o "QUAST" \\
        $args
        
    cp QUAST/report.tsv QUAST/report_rawassemblies.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}
