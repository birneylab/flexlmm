// match samples in the relatedness matrix, covariates, and phenotypes, and exclude
// samples where either of those is missing

process MATCH_SAMPLES {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(pheno), val(pheno_name), path(null_design_matrix), path(perm_group), path(model_frame)

    output:
    tuple val(meta), path("*.K.rds"), path("*.y.rds"), path("*.X.rds") , emit: model_terms
    tuple val(meta), path("*.model_frame.rds"), path("*.y.rds"), path("*.X.rds") , emit: model_terms
    tuple val(meta), path("*.perm_group.matched.rds")                  , emit: perm_group
    tuple val(meta), path("*.sample.id")                               , emit: sample_ids
    path "versions.yml"                                                , emit: versions

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

    X <- readRDS("${null_design_matrix}")
    y <- readRDS("${pheno}")[,"${pheno_name}"]
    perm_group <- readRDS("${perm_group}")
    covar_qcovar <- readRDS("${covar_qcovar}")

    X <- X[apply(!is.na(X), all, MARGIN = 1), , drop = FALSE]
    y <- y[!is.na(y)]
    perm_group <- perm_group[!is.na(perm_group)]

    samples_X <- rownames(X)
    samples_y <- names(y)
    samples_perm_group <- names(perm_group)
    samples_covar_qcovar <- rownames(covar_qcovar)

    samples <- intersect(
        intersect(
            intersect(
                intersect(samples_K, samples_y),
                samples_X
            ), samples_perm_group
        ), samples_covar_qcovar
    )
    X <- X[match(samples, rownames(X)), , drop = FALSE]
    y <- y[match(samples, names(y))]
    K <- K[match(samples, rownames(K)), match(samples, colnames(K))]
    covar_qcovar <- covar_qcovar[match(samples, rownames(covar_qcovar)),]
    perm_group <- perm_group[match(samples, names(perm_group))]

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(X)))
    stopifnot(all(names(y) == rownames(K)))
    stopifnot(all(names(y) == colnames(K)))
    stopifnot(all(names(y) == names(perm_group)))
    stopifnot(all(names(y) == rownames(covar_qcovar)))
    stopifnot(
        sum(is.na(K)) + sum(is.na(C)) + sum(is.na(y)) + sum(is.na(perm_group)) + sum(is.na(covar_qcovar)) == 0
    )

    message(length(samples), " sample intersect in all sets and have no missing values")

    write.table(
        samples,
        "${prefix}.sample.id",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
    )
    saveRDS(K, "${prefix}.K.rds")
    saveRDS(y, "${prefix}.y.rds")
    saveRDS(X, "${prefix}.X.rds")
    saveRDS(perm_group, "${prefix}.perm_group.matched.rds")

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
    touch ${prefix}.K.rds
    touch ${prefix}.y.rds
    touch ${prefix}.X.rds
    touch ${prefix}.sample.id
    touch ${prefix}.perm_group.matched.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
