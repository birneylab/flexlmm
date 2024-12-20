process TRANSFORM_PHENOTYPES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam), path(pheno)

    output:
    tuple val(meta), path("*.cov") , emit: pheno
    path "versions.yml"            , emit: versions

    when:
    params.quantile_normalise || params.standardise

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
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --pheno $pheno \\
        --pgen $pgen \\
        --psam $psam \\
        --pvar $pvar \\
        $cmd \\
        --write-covar 'cols=maybefid,maybesid,phenos' \\
        --out $prefix \\
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
    touch ${prefix}.cov

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
