// match samples in the relatedness matrix, covariates, and phenotypes, and exclude
// samples where either of those is missing

process MATCH_SAMPLES {
    tag "$meta.id"
    label 'process_low'

    // mulled r-data.table
    conda "bioconda::mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68==f20c3bc5c88913df9b835378643ab86f517a3dcf-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' :
        'biocontainers/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' }"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(pheno), val(pheno_name), path(null_design_matrix)

    output:
    tuple val(meta), path("*.K.rds"), path("*.y.rds"), path("*.C.rds") , emit: model_terms
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

    C <- readRDS("${null_design_matrix}")
    y <- readRDS("${pheno}")[,"${pheno_name}"]

    C <- C[apply(!is.na(C), all, MARGIN = 1),]
    y <- y[!is.na(y)]

    samples_C <- rownames(C)
    samples_y <- names(y)

    samples <- intersect(intersect(samples_K, samples_y), samples_C)
    C <- C[match(samples, rownames(C)),]
    y <- y[match(samples, names(y))]

    stopifnot(names(y) == rownames(C))
    stopifnot(names(y) == rownames(K))
    stopifnot(names(y) == colnames(K))
    stopifnot(sum(is.na(K)) + sum(is.na(C)) + sum(is.na(y)) == 0)

    message(length(samples), "sample intersect in all sets and have no missing values")

    saveRDS(K, "${prefix}.K.rds")
    saveRDS(y, "${prefix}.y.rds")
    saveRDS(C, "${prefix}.C.rds")

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
    touch ${prefix}.C.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
