// From each permutation collects the minimum p-value and returns a vector of such
// p-values with length equal to the number of permutations
process GET_MIN_P_DISTRIBUTION {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(perms)
    val nperms

    output:
    tuple val(meta), path("*.min_p_dist.rds") , emit: min_p_dist
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    p_vec  <- numeric(0)
    n_snps_vec <- numeric(0)
    perm_vec <- numeric(0)

    if (${nperms} > 1) pb <- txtProgressBar(1, ${nperms}, style = 3)

    for ( i in 1:${nperms} ) {
        if (${nperms} > 1) setTxtProgressBar(pb, i)
        all_files <- list.files(pattern = sprintf("perm%s.tsv.gwas.gz", i))
        min_p <- 2
        n_snps <- 0
        for ( perm_file in all_files ) {
            con <- gzfile(perm_file, "r")
            header <- readLines(con, n = 1)
            colnames <- strsplit(header, "\\t")[[1]]
            p_index <- which(colnames == "pval")
            stopifnot(length(p_index) == 1)

            line <- readLines(con, n = 1)
            while ( length(line) > 0 ) {
                new_p <- as.numeric(strsplit(line, "\\t")[[1]][[p_index]])
                min_p <- min(min_p, new_p, na.rm = TRUE)
                n_snps <- n_snps + 1
                line <- readLines(con, n = 1)
            }
            close(con)
        }
        p_vec <- c(p_vec, min_p)
        n_snps_vec <- c(n_snps_vec, n_snps)
        perm_vec <- c(perm_vec, i)
    }

    if (${nperms} > 1) close(pb)

    stopifnot(length(p_vec) == ${nperms})
    stopifnot(perm_vec == 1:${nperms})
    stopifnot(length(n_snps_vec) == ${nperms})
    stopifnot(all(p_vec <= 1 | p_vec == 2))
    stopifnot(all(p_vec >= 0))
    stopifnot(all(!is.na(n_snps_vec)))
    stopifnot(length(unique(n_snps_vec)) == 1)

    df <- data.frame(permutation = perm_vec, min_p = p_vec, n_snps = n_snps_vec)
    saveRDS(df, "${prefix}.min_p_dist.rds")

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
    touch ${prefix}.min_p_dist.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
