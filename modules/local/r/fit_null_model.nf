process FIT_NULL_MODEL {
    tag "$meta.id"
    label 'process_low'

    container 'saulpierotti-ebi/pgenlibr@sha256:0a606298c94eae8d5f6baa76aa1234fa5e7072513615d092f169029eacee5b60'

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
    y <- l[["y"]]
    X <- l[["X"]]

    fit <- .lm.fit(x = X, y = y)
    ll <- stats:::logLik.lm(fit)
    beta <- stats::coef(model)
    e <- stats::resid(fit)

    # the error is uncorrelated since I regressed the covariance structure but not
    # the realised residuals, which have correlation structure V = I - X * (t(X) * X)^(-1) * t(X)
    # See Abney, Genetic Epidemiology, Vol. 39, No. 4, 249–258, 2015.

    # this is the resiudal covariance matrix V
    I <- diag(length(y))
    V <- I - X %*% solve(t(X) %*% X, t(X))

    # V is symmetric and idempotent (V^2 = V) because the error covariance structure has been
    # regressed out already, so all the eigenvalues are either 0 or 1. 
    # do SVD and keep only eigenvalues of 1
    # since eigenvalues are 1 sqrt(V) = V
    V.eig <- eigen(V)
    V <- V.eig[["vectors"]][, V.eig[["values"]] > 0.9]
    
    # decorrelate the residual vector to get permutable residuals
    e.p <- crossprod(V, e)

    saveRDS(list(e = e.p, ll = ll, y.pred = fitted(fit)), "${prefix}.null_model.rds")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_pgenlibr <- utils::packageVersion("pgenlibr")
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    pgenlibr: %s", ver_pgenlibr),
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
        pgenlibr: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"pgenlibr\\")))")
    END_VERSIONS
    """
}
