process GREML {
    // gcta needs FID as first column, absent from pheno and qcovar, and no header
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gcta=1.94.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta) , path(grm), path(grm_id), val(pheno_idx)
    tuple val(meta2), path(pheno)
    tuple val(meta3), path(qcovar)

    output:
    tuple val(meta), path("*.hsq") , emit: hsq
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def mem_mb        = task.memory.toMega()
    def create_qcovar = qcovar ?
        "cat ${qcovar} | tail -n +2 | sed 's/^/0\\t/' > ${qcovar}.reformatted" :
        ""
    def qcovar_cmd = qcovar ? "--qcovar ${qcovar}.reformatted" : ""
    """
    cat $pheno | tail -n +2 | sed 's/^/0\\t/' > ${pheno}.reformatted
    ${create_qcovar}

    gcta64 \\
        --threads $task.cpus \\
        --grm-bin $grm.simpleName \\
        --pheno ${pheno}.reformatted \\
        --mpheno ${pheno_idx+1} \\
        $qcovar_cmd \\
        --out $prefix \\
        --reml \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtca: \$(gcta64 | grep version | sed -e 's/.*v//' -e 's/ .*//')
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
        gtca: \$(gcta64 | grep version | sed -e 's/.*v//' -e 's/ .*//')
    END_VERSIONS
    """
}
