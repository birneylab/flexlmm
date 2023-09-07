process GET_DESIGN_MATRIX {
    // Both covar and qcovar are exported to GCTA qcovars because factors are converted
    // to quantitative 0-1 dummy encoded variables
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::gcta=1.94.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gcta:1.94.1--h9ee0642_0' :
        'biocontainers/gcta:1.94.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(covar)
    tuple val(meta), path(qcovar)
    val formula_raw

    output:
    tuple val(meta), path("*.pheno" ) , emit: pheno
    tuple val(meta), path("*.covar" ) , emit: covar
    tuple val(meta), path("*.qcovar") , emit: qcovar
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def formula =
    """
    #!/usr/bin/env Rscript

    covar <- read.table(
        "${covar}" ,
        header = TRUE,
        comment.char = "",
        check.names = FALSE,
        colClasses = "factor",
    )

    qcovar_names <- names(read.table(
        "${qcovar}",
        header = TRUE,
        comment.char = "",
        check.names = FALSE,
        nrows = 1
    ))

    qcovar_classes <- sapply(
        qcovar_names,
        function(x){if (x == "#ID") "factor" else "numeric"}
    )

    qcovar <- read.table(
        "${qcovar}",
        header = FALSE,
        col.names = qcovar_names,
        colClasses = qcovar_classes
    )

    df <- merge(covar, qcovar, by = "#ID")
    C <- model.matrix(${formula}, data = df)
    out <- data.frame(ID = df["#ID"], FID = 0, C)
    write.table(out, "${prefix}.gcta_qcovar.tsv", sep = "\\t")
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def mem_mb      = task.memory.toMega()
    """
    touch ${prefix}.gcta_qcovar.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
