process MANHATTAN {
    tag "$meta.id"
    label 'process_low'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(gwas), path(min_p)
    val p_thr

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plot
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "png"
    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    gwas_files <- list.files(pattern = "*.gwas.tsv.gz")
    df <- lapply(gwas_files, read_tsv) %>% bind_rows()
    df_p_min <- readRDS("${min_p}")
    perm_thr <- as.numeric(quantile(df_p_min[["min_p"]], ${p_thr}))
    bonferroni_thr <- ${p_thr} / nrow(df)

    cumsums <- df %>%
        group_by(chr) %>%
        summarise(max_pos = max(pos)) %>%
        mutate(pos_add = lag(cumsum(max_pos), default = 0)) %>%
        select(chr, pos_add)

    df <- df %>%
        inner_join(cumsums, by = "chr") %>%
        mutate(pos_cum = pos + pos_add)

    axis_set <- df %>%
        group_by(chr) %>%
        summarize(center = mean(pos_cum))

    p <- ggplot(df, aes(x = pos_cum, y = -log10(lrt_p), color = as.factor(chr))) +
        geom_point(size = 0.8, alpha = 0.75) +
        geom_hline(yintercept = -log10(perm_thr), color = "red") +
        geom_hline(yintercept = -log10(bonferroni_thr), color = "blue") +
        theme_minimal_hgrid(18) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_x_continuous(label = axis_set[["chr"]], breaks = axis_set[["center"]]) +
        scale_color_manual(
            values = rep(c("black", "gray"), unique(length(axis_set[["chr"]])))
        ) +
        labs(x = NULL, y = bquote("-log"[10]~italic(p))) +
        theme(legend.position = "none") +
        annotate(
            geom = "text",
            y = -log10(bonferroni_thr),
            x = Inf,
            label = "Bonferroni ${p_thr}",
            size = 6,
            vjust = ifelse(-log10(bonferroni_thr) >= -log10(perm_thr), -0.2, 1.2),
            hjust = 1
        ) +
        annotate(
            geom = "text",
            y = -log10(perm_thr),
            x = Inf,
            label = "Permutations ${p_thr}",
            size = 6,
            vjust = ifelse(-log10(perm_thr) > -log10(bonferroni_thr), -0.2, 1.2),
            hjust = 1
        ) +
        ggtitle("${prefix}")

    ggsave("${prefix}.${suffix}", width = 30, height = 6, bg = "white")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_tidyverse <- utils::packageVersion("tidyverse")
    ver_cowplot <- utils::packageVersion("cowplot")

    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-tidyverse: %s", ver_tidyverse),
            sprintf("    r-cowplot: %s", ver_cowplot),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-tidyverse: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"tidyverse\\")))")
        r-cowplot: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"cowplot\\")))")
    END_VERSIONS
    """
}

process QQ {
    tag "$meta.id"
    label 'process_low'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(gwas)

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plot
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "png"
    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    gwas_files <- list.files(pattern = "*.gwas.tsv.gz")
    df <- lapply(gwas_files, read_tsv) %>%
        bind_rows() %>%
        arrange(lrt_p) %>%
        reframe(
            sample = -log10(lrt_p),
            theoretical = -log10(qunif(ppoints(lrt_p)))
        )

    p <- ggplot(df, aes(x = theoretical, y = sample)) +
        geom_point(size = 1) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme_cowplot(18) +
        labs(y = "Sample", x = "Theoretical")

    ggsave("${prefix}.${suffix}", p, bg = "white")

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_tidyverse <- utils::packageVersion("tidyverse")
    ver_cowplot <- utils::packageVersion("cowplot")

    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-tidyverse: %s", ver_tidyverse),
            sprintf("    r-cowplot: %s", ver_cowplot),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    touch ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-tidyverse: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"tidyverse\\")))")
        r-cowplot: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"cowplot\\")))")
    END_VERSIONS
    """
}


process RELATEDNESS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-complexheatmap==2.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-complexheatmap:2.16.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-complexheatmap:2.16.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(K), path(samples)

    output:
    tuple val(meta), path("${prefix}.png") , emit: plot
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library("ComplexHeatmap")
    library("circlize")

    samples <- read.table("${samples}")[[1]]
    K_vec <- readBin("${K}", what = "double", n = (length(samples))**2)
    K <- matrix(K_vec, nrow = length(samples))
    K.clust <- hclust(as.dist(max(K) - K))

    png("${prefix}.png", width = 9.4, height = 7)
    Heatmap(
        K,
        cluster_rows = K.clust,
        cluster_columns = K.clust,
        col = colorRamp2(c(min(K), mean(K), max(K)), c("blue", "white", "red")),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        heatmap_legend_param = list(title = "Relatedness"),
        use_raster = FALSE
    )
    dev.off()

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_complexheatmap <- utils::packageVersion("ComplexHeatmap")

    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-ComplexHeatmap: %s", ver_complexheatmap),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-ComplexHeatmap: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"ComplexHeatmap\\")))")
    END_VERSIONS
    """
}
