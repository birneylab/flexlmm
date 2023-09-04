process JOIN_RSQUARE_CHR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-data.table=1.14.8  r-r.utils=2.12.2"
    container "saulpierotti-ebi/r_datatable_tidyverse:0.1"

    input:
    tuple val(meta), path(rsquare_per_site)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: rsquare_per_site
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
    """
    #!/usr/bin/env Rscript

    library("data.table")

    setDTthreads(${task.cpus})

    l_names <- list.files(pattern = "*.txt.gz")
    l_df <- lapply(l_names, fread, header = TRUE, sep = '\t')
    df <- rbindlist(l_df)

    fwrite(df, "${prefix}.tsv.gz", sep = '\t')

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
        r-data.table: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"data.table\\")))")
    END_VERSIONS
    """
}
