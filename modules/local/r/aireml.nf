// determine variance components from the relatedness matrix, phenotypes, and covariates
process AIREML {
    // Both covar and qcovar are exported to GCTA qcovars because factors are converted
    // to quantitative 0-1 dummy encoded variables
    tag "$meta.id"
    label 'process_low'

    container "saulpierotti-ebi/gaston"

    input:
    tuple val(meta ), path(grm_bin), path(grm_id), path(pheno), val(pheno_name)
    tuple val(meta2), path(null_design_matrix)

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

    samples <- read.table("${grm_id}", header = FALSE, check.names = FALSE)[,1]
    K <- matrix(
        readBin("${grm_bin}", what="numeric", n=length(samples)**2),
        ncol = length(samples)
    )
    colnames(K) <- samples
    rownames(K) <- samples

    C <- readRDS("${null_design_matrix}")
    C <- C[match(samples, rownames(C)),]
    stopifnot(all(samples == rownames(C)))

    y <- readRDS("${pheno}")[,"${pheno_name}"]
    y <- y[match(samples, names(y))]
    stopifnot(all(samples == names(y)))

    fit <- gaston::lmm.aireml(y, C, K, verbose = TRUE)
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
    touch ${prefix}.covariate_mat.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-gaston: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"gaston\\")))")
    END_VERSIONS
    """
}
