process GET_DESIGN_MATRIX {
    // Both covar and qcovar are exported to GCTA qcovars because factors are converted
    // to quantitative 0-1 dummy encoded variables
    tag "$meta.id"
    label 'process_low'

    // mulled r-data.table
    conda "bioconda::mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68==f20c3bc5c88913df9b835378643ab86f517a3dcf-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' :
        'biocontainers/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' }"

    input:
    tuple val(meta), path(covar), path(qcovar)
    path pheno
    val null_model_formula

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

    library("data.table")
    setDTthreads(${task.cpus})

    pheno <- readRDS("${pheno}")
    clean_colnames <- function(n){gsub("#", "", n)}

    if ( ${load_covar} ){
        covar <- fread("${covar}" ,
            header = TRUE,
            sep = "\\t",
            check.names = FALSE,
            colClasses = "character"
        )

        colnames(covar) <- clean_colnames(colnames(covar))
        if ("FID" %in% colnames(covar)) covar[, FID := NULL]
    }

    if ( ${load_qcovar} ){
        qcovar <- fread(
            "${qcovar}",
            header = TRUE,
            sep = "\\t",
            check.names = FALSE,
            colClasses = "character"
        )

        colnames(qcovar) <- clean_colnames(colnames(qcovar))
        if ("FID" %in% colnames(qcovar)) qcovar[, FID := NULL]

        qcovar <- cbind(
            qcovar[, .(IID)],
            qcovar[, lapply(.SD, as.numeric), .SDcols = !"IID"]
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
        df <- data.table(IID = pheno[,"IID"])
    }

    null_model <- formula(${null_model_formula})
    null_model_RHS <- update(null_model, NULL ~ .)
    C <- model.matrix(null_model_RHS, data = df)
    rownames(C) <- df[["IID"]]
    saveRDS(C, "${prefix}.covariate_mat.rds")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_datatable <- utils::packageVersion("data.table")
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-data.table: %s", ver_datatable),
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
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"data.table\\")))")
    END_VERSIONS
    """
}
