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
    path fixed_effects_formula
    val  permute_by

    output:
    tuple val(meta), path("*.covariate_mat.rds") , emit: C
    tuple val(meta), path("*.gxe_frame.rds")     , emit: gxe_frame
    tuple val(meta), path("*.perm_group.rds")    , emit: perm_group
    path "versions.yml"                          , emit: versions

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

    covariate_formula <- readRDS("${covariate_formula}")
    fixed_effects_formula <- readRDS("${fixed_effects_formula}")
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

    fixed_vars <- all.vars(fixed_effects_formula)
    # this includes any covariate used (even if not in a gxe term) but it is fine, it
    # will just not be used downstream
    gxe_vars <- fixed_vars[fixed_vars != "x"]
    if ( length(gxe_vars) > 0 ){
        gxe_frame <- subset(df, select = gxe_vars)
        rownames(gxe_frame) <- df[["IID"]]
    } else {
        gxe_frame <- data.frame(row.names = samples)
    }
    saveRDS(gxe_frame, "${prefix}.gxe_frame.rds")

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
    saveRDS(perm_group, "${prefix}.perm_group.rds")

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
    touch ${prefix}.gxe_frame.rds
    touch ${prefix}.perm_group.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
