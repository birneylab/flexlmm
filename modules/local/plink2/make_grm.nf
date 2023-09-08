process MAKE_GRM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta) , path(pgen), path(psam), path(pvar), val(chr_exclude)
    tuple val(meta2), path(freq)

    output:
    tuple val(meta), path("*.grm.bin"), path("*.grm.id"), path("*.grm.N.bin") , emit: grm
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args   ?: ''
    def args2       = task.ext.args2  ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def mem_mb      = task.memory.toMega()
    def exclude_cmd = chr_exclude ? "--not-chr ${chr_exclude}" : ""
    def freq_cmd    = freq ? "--read-freq ${freq}" : ""
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --pgen $pgen \\
        --psam $psam \\
        --pvar $pvar \\
        ${exclude_cmd} \\
        ${freq_cmd} \\
        --out $prefix \\
        ${args} \\
        --make-grm-bin $args2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args   ?: ''
    def args2       = task.ext.args2  ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def mem_mb      = task.memory.toMega()
    def exclude_cmd = chr_exclude ? "--not-chr ${chr_exclude}" : ""
    def freq_cmd    = freq ? "--read-freq ${freq}" : ""
    """
    touch ${prefix}.grm.bin
    touch ${prefix}.grm.id
    touch ${prefix}.grm.N.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
