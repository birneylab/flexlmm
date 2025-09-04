process FIT_NULL_MODEL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(mm)

    output:
    tuple val(meta), path("*.null_model.rds") , emit: null_model
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    l <- readRDS("${mm}")
    y.mm <- l[["y.mm"]]
    X.mm <- l[["X.mm"]]

    stopifnot(all(!is.null(names(y.mm))))
    stopifnot(all(names(y.mm) == rownames(X.mm)))
    stopifnot(sum(is.na(X.mm)) + sum(is.na(y.mm)) == 0)

    fit <- lm.fit(x = X.mm, y = y.mm)
    ll <- stats:::logLik.lm(fit)
    e <- stats::resid(fit)

    # the error term is uncorrelated since I regressed the covariance structure in the DECORRELATE
    # process but not the realised residuals, which have correlation structure
    #
    # V = I - X * (t(X) * X)^(-1) * t(X)
    #
    # V is also called annihilator matrix and H = X * (t(X) * X)^(-1) * t(X) is the hat or projection matrix
    # that translates observed y to predicted y (called y_hat)
    #
    # See Abney, Genetic Epidemiology, Vol. 39, No. 4, 249–258, 2015.

    # this is the resiudal covariance matrix V
    I <- diag(length(y.mm))
    H <- X.mm %*% solve(t(X.mm) %*% X.mm, t(X.mm)) # hat matrix
    V <- I - H

    # V is symmetric and idempotent (V^2 = V) so all the eigenvalues are either 0 or 1. 
    # I do SVD and keep only eigenvalues of 1 and since eigenvalues are 1 sqrt(V) = V
    # U1 is V with the 0 eigenvalues and corresponding eigenvectors removed
    # if symmetric = FALSE, numerical errors can cause small complex eigenvalues
    # and fail when filtering at > 0.9
    V.eig <- eigen(V, symmetric = TRUE)
    U1 <- V.eig[["vectors"]][, V.eig[["values"]] > 0.9]
    
    # decorrelate the residual vector to get permutable residuals
    e.p <- crossprod(U1, e)

    # I save the elements needed to permute the phenotypes later on:
    # e.p, y.pred, and U1 because after permuting e.p I can get a new y as
    # y = y.pred + U1 %*% e.p
    #
    # I also save the null model log likelyhood to compute p-values later on
    saveRDS(list(e.p = e.p, y.mm.pred = fitted(fit), U1 = U1, ll = ll), "${prefix}.null_model.rds")

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
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.null_model.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
