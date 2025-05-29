process GET_VAR_IDX_EQTL {
    tag "$meta.id"
    label 'process_medium'

    conda "r-base=4.3.1 r-tidyverse=2.0.0"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam), val(chr)
    val(pheno_names)
    val(window)
    path(gtf)


    output:
    tuple val(meta), path("*.var_idx.rds") , emit: var_idx
    path "versions.yml"                    , emit: versions
    tuple val(meta), path("*_chr_pheno_map.txt"), emit: chr_pheno_map

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library("tidyverse")

    gtf <- read.delim("${gtf}", header = FALSE, comment.char = "#", sep = "\t")
    colnames(gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attributes")

    gene_names <- c(${pheno_names.collect { "\"${it}\"" }.join(", ")})

    genes <- gtf %>%
        filter(type == "gene" & chr == "${chr}") %>%
        mutate(
           gene_id = str_extract(attributes, "gene_id [^;]+") %>%
           str_replace("gene_id ", "")
        ) %>%
        filter(gene_id %in% gene_names)
    print(genes)

    chr_pheno_map <- genes %>%
        select(chr, gene_id) %>%
        distinct()

    cat("chr_pheno_map:\\n")
    print(chr_pheno_map)

    write.table(
        chr_pheno_map,
        file = "${prefix}_chr_pheno_map.txt",
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )

    pvar <- read.table(
        text = sub("^#", "", readLines("${pvar}")), header = TRUE
    )

    # identify SNPs within the gene window
    snps_in_window <- lapply(seq_len(nrow(genes)), function(i) {
        gene <- genes[i, ]
        snps <- which(
            pvar[["CHROM"]] == "${chr}" &
            pvar[["POS"]] >= (gene[["start"]] - as.numeric(${window})) &
            pvar[["POS"]] <= (gene[["end"]] + as.numeric(${window}))
        )
        # if no SNPs are found, return an empty data frame with the gene_id
       if (length(snps) == 0) {
        return(data.frame(gene_id = as.character(gene[["gene_id"]]), snps = NA_integer_))
       }
        # if SNPs are found, return the SNPs data frame
        data.frame(gene_id = as.character(gene[["gene_id"]]), snps = snps)  
    }) %>%
    bind_rows()

    saveRDS(snps_in_window, "${prefix}.var_idx.rds")

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
}

