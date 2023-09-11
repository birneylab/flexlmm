process PHENO_TO_RDS {
    tag "$meta.id"
    label 'process_low'

    // mulled r-data.table
    conda "bioconda::mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68==f20c3bc5c88913df9b835378643ab86f517a3dcf-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' :
        'biocontainers/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' }"

    input:
    tuple val(meta), path(pheno)

    output:
    tuple val(meta), path("*.pheno.rds")       , emit: pheno
    tuple val(meta), path("*.pheno_names.txt") , emit: pheno_names
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

    pheno <- fread(
        "${pheno}" ,
        header = TRUE,
        sep = "\\t",
        check.names = FALSE,
        colClasses = "character"
    )

    clean_colnames <- function(n){gsub("#", "", n)}
    colnames(pheno) <- clean_colnames(colnames(pheno))

    if ("FID" %in% colnames(pheno)) pheno[, FID := NULL]

    Y <- as.matrix(pheno[, lapply(.SD, as.numeric), .SDcols = !"IID"])
    rownames(Y) <- pheno[["IID"]]

    fwrite(
        data.table(colnames(Y)),
        "${prefix}.pheno_names.txt",
        col.names = FALSE
    )
    saveRDS(Y, "${prefix}.pheno.rds")

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
    touch ${prefix}.pheno.rds
    echo "stub_pheno" > ${prefix}.pheno_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"data.table\\")))")
    END_VERSIONS
    """
}
