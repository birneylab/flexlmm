process MANHATTAN {
    tag "$meta.id"
    label 'process_low'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(gwas), path(min_p)
    val p_thr

    output:
    tuple val(meta), path("*.pdf") , emit: plot
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
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
        geom_point(size = 1, alpha = 0.75) +
        geom_hline(yintercept = -log10(perm_thr), color = "red") +
        geom_hline(yintercept = -log10(bonferroni_thr), color = "blue") +
        theme_minimal_hgrid(18) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
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
            vjust = -1,
            hjust = 1
        ) +
        annotate(
            geom = "text",
            y = -log10(perm_thr),
            x = Inf,
            label = "Permutations ${p_thr}",
            size = 6,
            vjust = -1,
            hjust = 1
        ) +
        ggtitle("${prefix}")

    ggsave("${prefix}.${suffix}", p)

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
    tuple val(meta), path("*.pdf") , emit: plot
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    gwas_files <- list.files(pattern = "*.gwas.tsv.gz")
    df <- lapply(gwas_files, read_tsv) %>%
        bind_rows() %>%
        reframe(
            sample = -log10(lrt_p),
            theoretical = -log10(rank(lrt_p) / n())
        )

    p <- ggplot(df, aes(x = theoretical, y = sample)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme_cowplot(18) +
        labs(y = "Sample", x = "Theoretical")
        ggsave("${prefix}.${suffix}", p)

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
