process GET_CHR_NAMES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(pvar)

    output:
    tuple val(meta), path("*.txt") , emit: out
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat $pvar    | \\
    grep -v '^#' | \\
    cut -f1      | \\
    sort -u      > \\
    ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zstd: \$(zstd --version | sed 's/.*v\\(.*\\),.*/\\1/')
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def mem_mb      = task.memory.toMega()
    def exclude_cmd = chr_exclude ? "--not-chr ${chr_exclude}" : ""
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zstd: \$(zstd --version | sed 's/.*v\\(.*\\),.*/\\1/')
    END_VERSIONS
    """
}
