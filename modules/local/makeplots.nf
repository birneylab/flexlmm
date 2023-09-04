process PLOT_INFO_SCORE {
    tag "$meta.id"
    label 'process_single'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(info_score)

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plots
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    #! /usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    df <- read_csv("${info_score}")

    p <- ggplot(df, aes(x=as.numeric(info_score))) +
        labs(x = "Info Score", y = "Cumulative Density") +
        stat_ecdf(geom = "step") +
        theme_cowplot(18) +
        scale_x_continuous(limits = c(0,1))

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
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

process PLOT_R2_SITES {
    tag "$meta.id"
    label 'process_single'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(r2_sites)

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plots
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    #! /usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    df <- read_tsv("${r2_sites}")

    p <- ggplot(df, aes(x=as.numeric(ds_r2))) +
        labs(x = bquote(Pearson~italic(r^2)), y = "Cumulative Density") +
        stat_ecdf(geom = "step") +
        theme_cowplot(18) +
        scale_x_continuous(limits = c(0,1))

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
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


process PLOT_R2_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(r2_samples)

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plots
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ""
    prefix   = task.ext.prefix ?: "${meta.id}"
    suffix   = task.ext.suffix ?: "pdf"
    """
    #! /usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    df <- read_delim("${r2_samples}", col_names = c("sample", "gt_r2", "ds_r2"))

    # see github.com/odelaneau/GLIMPSE/pull/180
    df["ds_r2"] = df["ds_r2"] ** 2

    p <- ggplot(df, aes(x = as.numeric(ds_r2), y = fct_reorder(sample, ds_r2))) +
        labs(x = bquote(Pearson~italic(r^2)), y = "Sample") +
        geom_col() +
        scale_x_continuous(limits = c(0,1)) +
        scale_y_discrete(
            labels = as_labeller( function(x){ substr(x, 1, 40) } )
        ) +
        theme_cowplot(18) +
        theme(axis.text.x = element_text(size = 12))

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
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


process PLOT_R2_MAF {
    tag "$meta.id"
    label 'process_single'

    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(r2_maf)

    output:
    tuple val(meta), path("${prefix}.${suffix}") , emit: plots
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ""
    prefix            = task.ext.prefix ?: "${meta.id}"
    suffix            = task.ext.suffix ?: "pdf"
    """
    #! /usr/bin/env Rscript

    library("tidyverse")
    library("cowplot")

    df <- read_delim(
        "${r2_maf}",
        col_names = c(
            "row_id",
            "n_genotypes",
            "mean_AF",
            "gt_r2",
            "ds_r2"
        )
    )

    p <- ggplot(df, aes(x = as.numeric(mean_AF), y = as.numeric(ds_r2))) +
        labs(y = bquote(Pearson~italic(r^2)), x = "Minor Allele Frequency") +
        geom_line() +
        geom_point(size = 3) +
        scale_y_continuous(limits = c(0,1)) +
        theme_cowplot(18)

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args   ?: ""
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
