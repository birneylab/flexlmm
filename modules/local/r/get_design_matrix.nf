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
    path model
    path null_model
    val  permute_by

    output:
    tuple val(meta), path("X_null.rds")      , emit: x_null
    tuple val(meta), path("perm_group.rds")  , emit: perm_group
    tuple val(meta), path("model_frame.rds") , emit: model_frame
    path "versions.yml"                      , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def prefix         = task.ext.prefix ?: "${meta.id}"
    def load_covar     = covar      ? "TRUE" : "FALSE"
    def load_qcovar    = qcovar     ? "TRUE" : "FALSE"
    def permute_by_set = permute_by ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript

    # [-2] extract the LHS only from the formula
    null_model <- readRDS("${null_model}")[-2]
    model <- readRDS("${model}")[-2]
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

        covar <- cbind(
            subset(covar, select = (colnames(covar) == "IID")),
            lapply(remove_iid(covar), as.factor)
        )
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
    rownames(df) <- df[["IID"]]
    df[["x"]] <- NA
    model_frame <- model.frame(model, data = df)
    saveRDS(model_frame, "model_frame.rds")

    if ( ${permute_by_set} ){
        permute_by <- '${permute_by}'
        if ( !(${load_covar}) ) stop("'permute_by' is set but no covariate file was provided")
        if ( !(permute_by %in% colnames(covar)) ) {
            stop("The value of 'permute_by' must be a column of the covariate file.")
        }
        perm_group <- df[,permute_by]
    } else {
        perm_group <- rep(1, nrow(df))
    }
    names(perm_group) <- df[["IID"]]
    saveRDS(perm_group, "perm_group.rds")

    X_null <- model.matrix(null_model, data = df)
    saveRDS(X_null, "X_null.rds")

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
    touch X_null.rds
    touch perm_group.rds
    touch model_frame.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
