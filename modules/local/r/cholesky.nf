// Performs the cholesky decomposition V=LL^t of the phenotype variance-covariance matrix
// V. V is calculated from the variance components estimates as:
//
// V = s2g K + s2e I
//
// Where K is the genetic relatedness matrix, s2g and s2e the variance components, and I
// the identity matrix.

process CHOLESKY {
    tag "$meta.id"
    label 'process_low'

    // mulled r-data.table
    conda "bioconda::mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68==f20c3bc5c88913df9b835378643ab86f517a3dcf-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' :
        'biocontainers/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' }"

    input:
    tuple val(meta), path(K), path(y), path(C), path(gaston_hsq)

    output:
    tuple val(meta), path("*.chol_L.rds") , emit: chol_L
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    K <- readRDS("${K}")
    y <- readRDS("${y}")
    C <- readRDS("${C}")

    stopifnot(names(y) == rownames(C))
    stopifnot(names(y) == rownames(K))
    stopifnot(names(y) == colnames(K))
    stopifnot(sum(is.na(K)) + sum(is.na(C)) + sum(is.na(y)) == 0)

    # Load variance components
    gaston_hsq <- readRDS("${gaston_hsq}")
    s2g <- gaston_hsq[["tau"]]
    s2e <- gaston_hsq[["sigma2"]]

    # phenotype variance/covariance matrix
    V <- s2g * K + diag(s2e, dim(K))
    L <- t(chol(V)) # R returns the upper Cholesky triangle
    colnames(L) <- colnames(K)
    rownames(L) <- rownames(K)

    saveRDS(L, "${prefix}.chol_L.rds")

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
    touch ${prefix}.chol_L.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
