// rotates the covariate matrix C and the phenotype vector y to decorrelate their components by
// solving with respect to the pre-computed cholesky factor L of its covariance matrix
process DECORRELATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(y), path(C), path(chol_L)

    output:
    tuple val(meta), path("*.y_mm.rds"), path("*.C_mm.rds") , emit: mm_rotation
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    L <- readRDS("${chol_L}")
    y <- readRDS("${y}")
    C <- readRDS("${C}")

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(C)))
    stopifnot(all(names(y) == rownames(L)))
    stopifnot(all(names(y) == colnames(L)))
    stopifnot(sum(is.na(L)) + sum(is.na(C)) + sum(is.na(y)) == 0)

    y.mm <- forwardsolve(L, y)
    C.mm <- forwardsolve(L, C)

    names(y.mm) <- names(y)
    rownames(C.mm) <- rownames(C)
    colnames(C.mm) <- colnames(C)

    saveRDS(y.mm, "${prefix}.y_mm.rds")
    saveRDS(C.mm, "${prefix}.C_mm.rds")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            "END_VERSIONS\\n",
            sep = "\\n"
        )
    )
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.y_mm.rds
    touch ${prefix}.C_mm.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
