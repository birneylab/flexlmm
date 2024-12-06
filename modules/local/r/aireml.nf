// determine variance components from the relatedness matrix, phenotypes, and covariates
process AIREML {
    tag "$meta.id"
    label 'process_low'

    container "saulpierotti-ebi/gaston"

    input:
    tuple val(meta), path(K), path(y), path(X)

    output:
    tuple val(meta), path("*.hsq.rds") , emit: hsq
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    K <- readRDS("${K}")
    y <- readRDS("${y}")
    X <- readRDS("${X}")

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(X)))
    stopifnot(all(names(y) == rownames(K)))
    stopifnot(all(names(y) == colnames(K)))
    stopifnot(sum(is.na(K)) + sum(is.na(X)) + sum(is.na(y)) == 0)

    fit <- gaston::lmm.aireml(y, X, K, verbose = TRUE) 
    saveRDS(fit, "${prefix}.hsq.rds")

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
    touch ${prefix}.hsq.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-gaston: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"gaston\\")))")
    END_VERSIONS
    """
}
