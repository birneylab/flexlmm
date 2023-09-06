process ESTIMATE_FREQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("plink2.afreq.zst") , emit: freq
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --vcf $vcf \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        --set-missing-var-ids @_#_\\\\\$r_\\\\\$a \\
        --freq zs ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch plink2.afreq.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
