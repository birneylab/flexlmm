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
    val null_model_formula

    output:
    tuple val(meta), path("*.gcta_qcovar.tsv") , emit: mat
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library("data.table")

    setDTthreads(${task.cpus})

    covar <- fread(
        "${covar}" ,
        header = TRUE,
        sep = "\\t",
        check.names = FALSE,
        colClasses = "character"
    )

    qcovar <- fread(
        "${qcovar}",
        header = TRUE,
        sep = "\\t",
        check.names = FALSE,
        colClasses = "character"
    )

    qcovar <- cbind(
        qcovar[, .(`#IID`)],
        qcovar[, lapply(.SD, as.numeric), .SDcols = !"#IID"]
    )

    null_model <- formula(${null_model_formula})
    formula    <- update(null_model, NULL ~ .)

    message(paste("Design matrix formula:", deparse(formula)))

    df <- merge(covar, qcovar, by = "#IID")
    C <- model.matrix(formula, data = df)
    out <- data.table(`#IID` = df[["#IID"]], C, check.names = FALSE)
    fwrite(
        out,
        "${prefix}.gcta_qcovar.tsv",
        sep = "\\t",
        col.names = TRUE
    )

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
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gcta_qcovar.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"data.table\\")))")
    END_VERSIONS
    """
}
