// concatenate permutation and chromosome minimum p-values
process CAT_PERM {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"
    
    input:
    tuple val(meta), path(perms)

    output:
    tuple val(meta), path("*.all_perms.tsv.gz") , emit: min_p_dist
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    echo "permutation\tchr\tmin_pval\tnvars" | gzip > header.txt.gz
    cat header.txt.gz ${perms} > ${prefix}.all_perms.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.all_perms.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}

// summarise by chromosome to get genome-wide permutation minimum p-value distributions
process SUMMARISE_BY_CHR {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"
    
    input:
    tuple val(meta), path(perm)
    val nperms

    output:
    tuple val(meta), path("*.genome_wide_min_p.tsv.gz") , emit: min_p_dist
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    df <- read.table("${perm}", header = TRUE)
    chrs <- unique(df[["chr"]])
    perms <- unique(df[["permutations"]])

    out <- aggregate(min_pval ~ permutation, data = df, FUN = min)
    stopifnot(nrow(out) == ${nperms})
    out_con <- gzfile("${prefix}.genome_wide_min_p.tsv.gz", "w")
    write.table(out, out_con, sep = "\\t", row.names = FALSE, col.names = TRUE)
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
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genome_wide_min_p.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
