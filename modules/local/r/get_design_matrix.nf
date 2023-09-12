// the design matrix for the null model, including the intercept if required
process GET_DESIGN_MATRIX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(covar), path(qcovar)
    path pheno
    path covariate_formula

    output:
    tuple val(meta), path("*.covariate_mat.rds") , emit: mat
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def load_covar  = covar  ? "TRUE" : "FALSE"
    def load_qcovar = qcovar ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript

    covariate_formula <- readRDS("${covariate_formula}")
    samples <- rownames(readRDS("${pheno}"))
    clean_colnames <- function(n){gsub("#", "", n)}
    remove_fid <- function(df){subset(df, select = (colnames(df) != "FID"))}
    remove_iid <- function(df){subset(df, select = (colnames(df) != "IID"))}

    if ( ${load_covar} ){
        covar <- read.table("${covar}" ,
            header = TRUE,
            sep = "\\t",
            check.names = FALSE,
            colClasses = "character",
            comment.char = ""
        )

        colnames(covar) <- clean_colnames(colnames(covar))
        covar <- remove_fid(covar)
    }

    if ( ${load_qcovar} ){
        qcovar <- read.table(
            "${qcovar}",
            header = TRUE,
            sep = "\\t",
            check.names = FALSE,
            colClasses = "character",
            comment.char = ""
        )

        colnames(qcovar) <- clean_colnames(colnames(qcovar))
        qcovar <- remove_fid(qcovar)

        qcovar <- cbind(
            subset(qcovar, select = (colnames(qcovar) == "IID")),
            lapply(remove_iid(qcovar), as.numeric)
        )
    }

    if ( ${load_covar} & ${load_qcovar} ){
        df <- merge(covar, qcovar, by = "IID")
    } else if ( ${load_covar} ){
        df <- covar
    } else if ( ${load_qcovar} ){
        df <- qcovar
    } else {
        # this is needed to still have an intercept with the right number of samples
        # in case of no covar and no qcovar
        df <- data.frame(IID = samples)
    }

    C <- model.matrix(covariate_formula, data = df)
    rownames(C) <- df[["IID"]]
    saveRDS(C, "${prefix}.covariate_mat.rds")

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
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.covariate_mat.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
