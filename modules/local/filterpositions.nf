// filter imputation positions accotding to previous imputation performance
process FILTER_POSITIONS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-data.table=1.14.8  r-r.utils=2.12.2"
    container "saulpierotti-ebi/r_datatable_tidyverse:0.1"

    input:
    tuple val(meta), path(performance_csv), val(performance_threshold)
    val filter_var

    output:
    tuple val(meta), path("*.tsv") , emit: posfile
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args   ?: ""
    """
    #!/usr/bin/env Rscript

    library("data.table")

    setDTthreads(${task.cpus})

    colnames <- c("chr", "pos", "ref", "alt")
    df <- fread("${performance_csv}", header = TRUE)
    df <- df[${filter_var} > ${performance_threshold}]
    df <- df[,..colnames]

    fwrite(df, "${prefix}.tsv", sep = "\\t", col.names = FALSE)

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_datatable <- utils::packageVersion("data.table")

    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-data.table: %s", ver_datatable),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\"data.table\")))")
    END_VERSIONS
    """
}


process REFORMAT_R2 {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-data.table=1.14.8  r-r.utils=2.12.2"
    container "saulpierotti-ebi/r_datatable_tidyverse:0.1"

    input:
    tuple val(meta), path(r2_per_site)

    output:
    tuple val(meta), path("*.csv.gz") , emit: r2
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args   ?: ""
    """
    #!/usr/bin/env Rscript

    library("data.table")

    setDTthreads(${task.cpus})

    df <- fread("${r2_per_site}", sep = "\t", header = TRUE)
    df <- df[, .(chr, pos, ref = allele1, alt = allele2, r2 = ds_r2)]

    fwrite(df, "${prefix}.csv.gz")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_datatable <- utils::packageVersion("data.table")

    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-data.table: %s", ver_datatable),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """
    touch ${prefix}.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\"data.table\")))")
    END_VERSIONS
    """
}
