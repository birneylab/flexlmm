process CHOLESKY {
    tag "$meta.id"
    label 'process_low'

    // mulled r-data.table
    conda "bioconda::mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68==f20c3bc5c88913df9b835378643ab86f517a3dcf-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' :
        'biocontainers/mulled-v2-a0002b961f72ad8f575ed127549e478f81093b68:f20c3bc5c88913df9b835378643ab86f517a3dcf-0' }"

    input:
    tuple val(meta), path(grm_bin), path(grm_id), path(grm_n), path(gcta_hsq)

    output:
    tuple val(meta), path("*.chol_L.rds") , emit: chol_L
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    ####################################################################################
    # from https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM
    ####################################################################################

    ReadGRMBin=function(prefix, AllN=F, size=4){
      sum_i=function(i){
        return(sum(1:i))
      }
      BinFileName=paste(prefix,".grm.bin",sep="")
      NFileName=paste(prefix,".grm.N.bin",sep="")
      IDFileName=paste(prefix,".grm.id",sep="")
      id = read.table(IDFileName)
      n=dim(id)[1]
      BinFile=file(BinFileName, "rb");
      grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
      NFile=file(NFileName, "rb");
      if(AllN==T){
        N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
      }
      else N=readBin(NFile, n=1, what=numeric(0), size=size)
      i=sapply(1:n, sum_i)
      return(list(diag=grm[i], off=grm[-i], id=id, N=N))
    }

    ####################################################################################

    GRM <- ReadGRMBin("${grm_bin.simpleName}")

    # make the matrix from the diagonal and lower triangle
    samples <- GRM[["id"]][["V2"]]
    n_samples <- length(samples)
    K <- matrix(NA, ncol = n_samples, nrow = n_samples)
    diag(K) <- GRM[["diag"]]
    K[lower.tri(K)] <- GRM[["off"]]
    K[upper.tri(K)] <- t(K)[upper.tri(K)]

    # Load variance components
    hsq_table <- read.table("${gcta_hsq}", nrows=3, header = TRUE)
    s2g <- hsq_table[hsq_table["Source"] == "V(1)", "Variance"]
    s2e <- hsq_table[hsq_table["Source"] == "V(e)", "Variance"]

    # phenotype variance/covariance matrix
    V <- s2g * K + diag(s2e, n_samples)
    L <- t(chol(V)) # R returns the upper Cholesky triangle
    colnames(L) <- samples
    rownames(L) <- samples

    saveRDS(L, "${prefix}.chol_L.rds")

    ver_sum_i <- function(i){r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            "END_VERSIONS",
            sep = "\\n"
        )
    )
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chol_L.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
