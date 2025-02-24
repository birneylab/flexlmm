process GENE_EXPR_TO_PHENO {
    label 'process_low'

    conda "r-base=4.3.1 r-tidyverse=2.0.0"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    path(pheno)

    output:
    path("phenotypes.tsv"), emit: gene_pheno
    path("versions.yml"), emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library("tidyverse")

    # Read the input expression table
    expression_data <- read.table("${pheno}", header = TRUE, sep = "\\t", check.names = FALSE)

    # Remove the `gene_name` column and retain only `gene_id` and sample values
    expression_data <- expression_data %>%
        select(-gene_name)

    # Transpose the expression table so samples are rows and genes are columns
    phenotype_data <- expression_data %>%
        column_to_rownames("gene_id") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "#IID")


    # Write the final phenotype file
    write.table(
        phenotype_data,
        "phenotypes.tsv",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    # Capture versions of R and required libraries
    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_tidyverse <- utils::packageVersion("tidyverse")

    # Save versions to a YAML file
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    r-tidyverse: %s", ver_tidyverse),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    """
    touch phenotypes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        r-tidyverse: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"tidyverse\\")))")
    END_VERSIONS
    """
}

