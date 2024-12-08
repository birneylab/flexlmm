// determine variance components from the relatedness matrix, phenotypes, and covariates
process AIREML {
    tag "$meta.id"
    label 'process_low'

    container "saulpierotti-ebi/gaston"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(pheno), val(pheno_name)
    tuple val(meta2), path(x_null)

    output:
    tuple val(meta), path("*.aireml.rds") , emit: aireml
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    samples_K <- read.table("${grm_id}", header = FALSE, check.names = FALSE)[,1]
    K <- matrix(
        readBin("${grm_bin}", what="numeric", n=length(samples_K)**2),
        ncol = length(samples_K)
    )
    colnames(K) <- samples_K
    rownames(K) <- samples_K
    stopifnot(sum(is.na(K)) == 0)

    y <- readRDS("${pheno}")[, "${pheno_name}"]
    X <- readRDS("${x_null}")

    y <- y[!is.na(y)]

    samples <- intersect(samples_K, names(y))
    X <- X[match(samples, rownames(X)), , drop = FALSE]
    y <- y[match(samples, names(y))]
    K <- K[match(samples, rownames(K)), match(samples, colnames(K))]

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(X)))
    stopifnot(all(names(y) == rownames(K)))
    stopifnot(all(names(y) == colnames(K)))
    stopifnot(sum(is.na(K)) + sum(is.na(X)) + sum(is.na(y)) == 0)

    fit <- gaston::lmm.aireml(y, X, K, verbose = TRUE) 
    saveRDS(list(fit = fit, K = K, y = y, X = X), "${prefix}.aireml.rds")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_gaston <- utils::packageVersion("gaston")
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-gaston: %s", ver_gaston),
            "END_VERSIONS\\n",
            sep = "\\n"
        )
    )
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aireml.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-gaston: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"gaston\\")))")
    END_VERSIONS
    """
}
