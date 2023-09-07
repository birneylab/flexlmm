process GREML {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gcta=1.94.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(grm), val(pheno_idx)
    tuple val(meta), path(pheno)
    tuple val(meta), path(covar)
    tuple val(meta), path(qcovar)

    output:
    tuple val(meta), path("*.hsq") , emit: hsq
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    def cmd    = params.quantile_normalise ?
        "--quantile-normalize" :
        params.standardise ?
        "--variance-standardize" :
        ""
    """
    gcta64 \\
        --threads $task.cpus \\
        --grm-bin $grm \\
        --pheno $pheno \\
        $covar_cmd \\
        $qcovar_cmd \\
        --out $prefix \\
        --reml \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def mem_mb      = task.memory.toMega()
    """
    touch ${prefix}.hsq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
