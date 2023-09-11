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

    pheno <- read.table(
        "${pheno}" ,
        header = TRUE,
        sep = "\\t",
        check.names = FALSE,
        colClasses = "character",
        comment.char = ""
    )

    clean_colnames <- function(n){gsub("#", "", n)}
    remove_fid <- function(df){subset(df, select = (colnames(df) != "FID"))}
    remove_iid <- function(df){subset(df, select = (colnames(df) != "IID"))}
    colnames(pheno) <- clean_colnames(colnames(pheno))
    pheno <- remove_fid(pheno)

    Y <- as.matrix(as.data.frame(lapply(remove_iid(pheno), as.numeric)))
    rownames(Y) <- pheno[,"IID"]

    write.table(
        colnames(Y),
        "${prefix}.pheno_names.txt",
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE
    )
    saveRDS(Y, "${prefix}.pheno.rds")

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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pheno.rds
    echo "stub_pheno" > ${prefix}.pheno_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
