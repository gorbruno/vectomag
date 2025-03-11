process MAKE_MMSEQS_TAXONOMY_TABLE {

    conda "conda-forge::python=3.12.8 conda-forge::pandas=2.2.3 conda-forge::xlsxwriter=3.2.2 conda-forge::biopython=1.85"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' :
        'quay.io/biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' }"

    input:
    tuple val(meta), path('mmseqs2/*')
    path ('contigs/*')
    val  outname

    output:
    tuple val(meta), path("*.csv")       , emit: csv
    tuple val(meta), path("*.xlsx")      , optional: true, emit: excel
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args       ?: ''
    def pattern = task.ext.pattern ?: ''
    """
    make_mmseqs_taxonomy_table.py \\
        --pattern_sample $pattern \\
        --mmseqs2_dir ./mmseqs2 \\
        --contigs_dir ./contigs \\
        --output_file $outname \\
        $args

    echo $meta.mmseqs2_db_name
    echo $meta

    if [[ $outname != "merged" ]]; then
        # find . -name "variants_long_table.*" -exec sh -c "mv \$1 ${outname}.variants.\${1##*.}" rename {} \; TODO
        mv variants_long_table.csv ${outname}.variants.csv
        mv variants_long_table.xlsx ${outname}.variants.xlsx #may fail
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
