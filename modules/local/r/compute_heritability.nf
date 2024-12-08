process COMPUTE_HERITABILITY {
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    path aireml

    output:
    path "heritability.tsv.gz" , emit: table
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    outname <- "heritability.tsv.gz"
    out_con <- gzfile(outname, "w")
    header <- "phenotype\\tgenetic_variance\\tresidual_variance\\tgrm_heritability"
    writeLines(header, out_con)

    f_list <- list.files(pattern = "\\\\.aireml\\\\.rds\$")
    for (f in f_list) {
        pheno <- gsub("\\\\.aireml\\\\.rds\$", "", f)
        pheno <- gsub("^[^_]*_", "", pheno)
        l <- readRDS(f)
        fit <- l[["fit"]]
        s2e <- fit[["sigma2"]]
        s2g <- fit[["tau"]]
        h2 <- s2g / (s2g + s2e)
        lineout <- sprintf(
            "%s\\t%s\\t%s\\t%s",
            pheno,
            s2g,
            s2e,
            h2
        )
        writeLines(lineout, out_con)
    }
   
    close(out_con)

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
    """
    touch heritability.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
