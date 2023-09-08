process FIT_MODEL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-rrbgen=0.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-rrbgen:0.0.6--r43h4ac6f70_10' :
        'biocontainers/r-rrbgen:0.0.6--r43h4ac6f70_10' }"

    input:
    tuple val(meta), path(chol_L), path(bgen), val(samplefile)

    output:
    tuple val(meta), path("*.pheno.rds") , emit: pheno
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library("data.table")

    setDTthreads(${task.cpus})

    L <- readRDS("${chol_L}")
    y <- fread("${pheno}", sep = "\\t", header = TRUE)[["${pheno_col}"]]

    y.mm <- forwardsolve(L, y)

    saveRDS(y.mm, "${prefix}.pheno.rds")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_datatable <- utils::packageVersion("data.table")
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-data.table: %s", ver_datatable),
            "END_VERSIONS\\n",
            sep = "\\n"
        )
    )
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.null_design_matrix.rds
    touch ${prefix}.pheno.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"data.table\\")))")
    END_VERSIONS
    """
}
