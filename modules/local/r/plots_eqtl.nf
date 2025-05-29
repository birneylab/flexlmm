process MANHATTAN {

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    path(gwas)  // Now takes final_combined_gwas.tsv.gz as input
    val p_thr

    output:
    path "manhattan_plot.png" , emit: plot
    path "versions.yml"       , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    # Load final combined GWAS file
    df <- read_tsv("${gwas}")

    # Compute Bonferroni threshold
    bonferroni_thr <- ${p_thr} / nrow(df)

    # Compute cumulative position
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

    # Manhattan plot
    p <- ggplot(df, aes(x = pos_cum, y = -log10(pval), color = as.factor(chr))) +
        geom_point(size = 0.8, alpha = 0.75) +
        geom_hline(yintercept = -log10(bonferroni_thr), color = "blue", linetype = "dotted") +
        theme_minimal() +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        scale_x_continuous(labels = axis_set[["chr"]], breaks = axis_set[["center"]]) +
        scale_color_manual(values = rep(c("black", "gray"), unique(length(axis_set[["chr"]])))) +
        labs(x = "Chromosome Position", y = bquote("-log"[10]~italic(p))) +
        theme(legend.position = "none") +
        ggtitle("Genome-wide Manhattan Plot (eQTL)")

    # Save plot
    ggsave("manhattan_plot.png", width = 15, height = 6, bg = "white")

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
}



process QQ {

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    path(gwas) 

    output:
    path "qq_plot.png" , emit: plot
    path "versions.yml" , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    # Load final combined GWAS file
    df <- read_tsv("${gwas}") %>%
        arrange(pval) %>%
        mutate(
            sample = -log10(pval),
            theoretical = -log10(qunif(ppoints(pval)))
        )
    
    # Calculate lambda (genomic inflation factor)
    df\$chi_sq <- qchisq(df\$pval, df = 1, lower.tail = FALSE)
    lambda <- median(df\$chi_sq) / 0.456  # 0.456 is expected median for df=1

    # QQ plot
    p <- ggplot(df, aes(x = theoretical, y = sample)) +
        geom_point(size = 1, alpha = 0.8) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        theme_minimal() +
        labs(y = "Observed -log10(p)", x = "Expected -log10(p)") +
        ggtitle("Q-Q Plot") +
        annotate("text", x = max(df\$theoretical) * 0.7, y = max(df\$sample) * 0.9, 
                 label = sprintf("λ = %.3f", lambda), 
                 size = 6, color = "blue") 

    # Save plot
    ggsave("qq_plot.png", width = 8, height = 8, bg = "white")

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
}

