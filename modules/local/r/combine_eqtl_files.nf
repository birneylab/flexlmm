process COMBINE_EQTL_FILES {

    
    conda "r-base=4.3.1 r-tidyverse=2.0.0 r-cowplot=1.1.1"
    container "saulpierotti-ebi/r_plotting:0.1"

    input:
    val gwas_files

    output:
    path "final_combined_gwas.tsv.gz" , emit: eqtl

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    gwas_files <- strsplit("${gwas_files}", " ")[[1]]

    combined_list <- list()

    for (file in gwas_files) {
      file_name <- basename(file)
      parts <- strsplit(file_name, "_")[[1]]
      gene <- strsplit(parts[length(parts)], "[.]")[[1]][1]

      df <- fread(cmd = paste("zcat", file), data.table = FALSE)

      if (nrow(df) == 0) next
      df\$gene <- gene
      combined_list[[gene]] <- df
    }

    combined_df <- do.call(rbind, combined_list)

    fwrite(combined_df, "final_combined_gwas.tsv.gz", sep = "\\t", quote = FALSE, row.names = FALSE)
    """
}

