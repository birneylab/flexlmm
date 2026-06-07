// Performs the cholesky decomposition V=LL^t of the phenotype variance-covariance matrix
// V. V is calculated from the variance components estimates as:
//
// V = s2g K + s2e I
//
// Where K is the genetic relatedness matrix, s2g and s2e the variance components, and I
// the identity matrix.
//
// Then, rotates the covariate matrix X and the phenotype vector y to decorrelate their components by
// solving with respect to the pre-computed cholesky factor L of the covariance matrix

process DECORRELATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(aireml)

    output:
    tuple val(meta), path("*.mm.rds") , emit: mm
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    l <- readRDS("${aireml}")
    fit <- l[["fit"]]
    K <- l[["K"]]
    y <- l[["y"]]
    X <- l[["X"]]

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(X)))
    stopifnot(all(names(y) == rownames(K)))
    stopifnot(all(names(y) == colnames(K)))
    stopifnot(sum(is.na(K)) + sum(is.na(X)) + sum(is.na(y)) == 0)

    # Load variance components
    s2g <- fit[["tau"]]
    s2e <- fit[["sigma2"]]

    # Correct K if it has negative eigenvalues (numerical noise in GRM)
    eig_K <- eigen(K, symmetric = TRUE)
    n_negative <- sum(eig_K\$values < 0)
    if (n_negative > 0) {
        min_eig_K <- min(eig_K\$values)
        stopifnot(min_eig_K >= ${params.min_allowed_eig})
        message(sprintf("Gene %s: K has %d negative eigenvalues (min: %.6e), applying ridge correction",
                        "${meta.id}", n_negative, min_eig_K))
        K <- K + diag(abs(min_eig_K) + 1e-6, nrow(K))
    }

    V <- s2g * K + diag(s2e, nrow(K))
    L <- try(t(chol(V)), silent = TRUE)

    if (inherits(L, "try-error")) {
        V.eig <- eigen(V, symmetric = TRUE)
        message(sprintf("Gene %s: Cholesky failed after K correction (V min eigenvalue: %.6e), correcting V",
                        "${meta.id}", min(V.eig\$values)))
        V.eig\$values[V.eig\$values <= 0] <- ${params.eig_replacement}
        V <- V.eig\$vectors %*% diag(V.eig\$values) %*% t(V.eig\$vectors)
        L <- t(chol(V))
    }

    colnames(L) <- colnames(K)
    rownames(L) <- rownames(K)

    y.mm <- forwardsolve(L, y)
    X.mm <- forwardsolve(L, X)

    names(y.mm) <- names(y)
    rownames(X.mm) <- rownames(X)
    colnames(X.mm) <- colnames(X)

    saveRDS(list(y.mm = y.mm, X.mm = X.mm, L = L), "${prefix}.mm.rds")

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
    touch ${prefix}.mm.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
