process FIT_MODEL {
    tag "$meta.id"
    label 'process_low'

    container 'saulpierotti-ebi/pgenlibr@sha256:0a606298c94eae8d5f6baa76aa1234fa5e7072513615d092f169029eacee5b60'

    input:
    tuple val(meta), path(mm), path(var_idx), path(null_model_fit), val(perm_seed)
    tuple val(meta2), path(pgen), path(pvar), path(psam)
    path model
    tuple val(meta3), path(model_frame)
    val use_dosage

    output:
    tuple val(meta), path("*.${perm_seed ? 'perm.tsv.gz' : 'tsv.gwas.gz'}") , emit: out
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def do_permute = perm_seed ? "TRUE" : "FALSE"
    def the_perm_seed  = perm_seed ?: "NULL"
    def pgenlibr_read_func = use_dosage ? "Read" : "ReadHardcalls"
    """
    #!/usr/bin/env Rscript

    l.mm <- readRDS("${mm}")
    var_idx <- readRDS("${var_idx}")
    l.null <- readRDS("${null_model_fit}")
    
    pvar_table <- read.table(
        # header starts with # and comment line with ##,
        # by removing the first # I make the header visible
        text = sub("^#", "", readLines("${pvar}")), header = TRUE
    )[var_idx,]
    pvar <- pgenlibr::NewPvar("${pvar}")
    pgen <- pgenlibr::NewPgen("${pgen}", pvar = pvar)
    psam <- read.table(
        "${psam}",
        sep = "\\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    clean_colnames <- function(n){
        n[n == "#IID"] <- "IID"
        n[n == "#FID"] <- "FID"
    }
    colnames(psam) <- clean_colnames(colnames(psam))
    buf <- pgenlibr::Buf(pgen)

    model <- readRDS("${model}")
    model_frame_raw <- readRDS("${model_frame}")

    y.mm <- l.mm[["y.mm"]]
    X.mm.null <- l.mm[["X.mm"]]
    L <- l.mm[["L"]]

    ll.null <- l.null[["ll"]]
    y.mm.pred <- l.null[["y.mm.pred"]]
    e.p <- l.null[["e.p"]]
    U1 <- l.null[["U1"]]
    
    # samples in y.mm are already intersected with psam[["IID"]] in aireml step when they are matched with the GRM, so no risk of NA here
    samples <- names(y.mm)
    pgen_order <- match(samples, psam[["IID"]])
    psam <- psam[pgen_order,]
    model_frame_raw <- model_frame_raw[match(samples, rownames(model_frame_raw)),]

    stopifnot(all(!is.null(names(y.mm))))
    stopifnot(all(names(y.mm) == rownames(X.mm.null)))
    stopifnot(all(names(y.mm) == rownames(L)))
    stopifnot(all(names(y.mm) == colnames(L)))
    stopifnot(all(names(y.mm) == psam[["IID"]]))
    stopifnot(all(names(y.mm) == names(y.mm.pred)))
    stopifnot(all(names(y.mm) == rownames(model_frame_raw)))
    stopifnot(
        (
            sum(is.na(L)) +
            sum(is.na(X.mm.null)) +
            sum(is.na(y.mm)) +
            sum(is.na(psam[["IID"]])) +
            sum(is.na(y.mm.pred)) +
            sum(is.na(model_frame_raw))
        ) == 0
    )
        
    # generate synthetic phenotype by permuting the uncorrelated transformation
    # of the residuals and fit the null model
    # null model fit is not needed if not permuting because we can use the one
    # obtained in the FIT_NULL_MODEL process
    if (${do_permute}) {
        set.seed(${the_perm_seed})
        y.mm <- y.mm.pred + drop(U1 %*% sample(e.p))
        stopifnot(all(rownames(X.mm.null) == rownames(y.mm)))
        fit <- .lm.fit(x = X.mm.null, y = y.mm)
        ll.null <- stats:::logLik.lm(fit)

        # an impossibly large value to initiate the variable
        min_pval <- 2
        # counter to be filled in in the for loop
        nvars_tested <- 0

        outname <- "${prefix}.perm.tsv.gz"
        out_con <- gzfile(outname, "w")
    } else {
        # generate SNP-wise output only for non-permuted version
        outname <- "${prefix}.tsv.gwas.gz"
        out_con <- gzfile(outname, "w")
        header <- "chr\\tpos\\tid\\tref\\talt\\tlrt_chisq\\tlrt_df\\tpval\\tbeta"
        writeLines(header, out_con)
    }

    nvars <- length(var_idx)
    pb <- txtProgressBar(1, nvars, style = 3)
    for (i in seq_along(var_idx)) {
        setTxtProgressBar(pb, i)
        the_var_idx <- var_idx[[i]]
        pgenlibr::${pgenlibr_read_func}(pgen, buf, the_var_idx)
        model_frame_raw[["x"]] <- buf[pgen_order]
        # needed to recompute expressions such as I(x == 1)
        model_frame <- model.frame(model[-2], data = model_frame_raw)
        X <- model.matrix(model[-2], data = model_frame)
        X.mm <- forwardsolve(L, X)
        colnames(X.mm) <- colnames(X)
        rownames(X.mm) <- rownames(X)

        stopifnot(all(rownames(X.mm) == rownames(y.mm)))
        fit <- .lm.fit(x = X.mm, y = y.mm)
        ll <- stats:::logLik.lm(fit)
        lrt_df <- attributes(ll)[["df"]] - attributes(ll.null)[["df"]]
        lrt_chisq <- 2 * as.numeric(ll - ll.null)
        pval <- (
            if (lrt_df > 0) pchisq(lrt_chisq, df = lrt_df, lower.tail = FALSE) else NA
        )

        if (${do_permute}) {
            nvars_tested <- nvars_tested + 1
            min_pval <- min(pval, min_pval, na.rm = TRUE)
        } else {
            id <- pgenlibr::GetVariantId(pvar, the_var_idx)
            chr <- pvar_table[i, "CHROM"]
            pos <- pvar_table[i, "POS"]
            ref <- pvar_table[i, "REF"]
            alt <- pvar_table[i, "ALT"]
            beta <- paste(
                colnames(X.mm),
                coef(fit),
                sep = "~",
                collapse = ","
            )
        
            lineout <- sprintf(
                "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",
                chr,
                pos,
                id,
                ref,
                alt,
                lrt_chisq,
                lrt_df,
                pval,
                beta
            )
            writeLines(lineout, out_con)
        }
    }

    pgenlibr::ClosePgen(pgen)
    close(pb)

    if (${do_permute}) {
        stopifnot(min_pval >= 0 && min_pval <= 1)
        stopifnot(nvars_tested == nvars)
        lineout <- sprintf("${the_perm_seed}\\t${meta.chr}\\t%s\\t%s", min_pval, nvars)
        writeLines(lineout, out_con)
        close(out_con)
    } else {
        close(out_con)
        # to make sure that the output has been written properly
        gwas <- read.table(outname, header = TRUE, sep = "\t")
        stopifnot(nrow(gwas) == nvars)
        stopifnot(ncol(gwas) == 9)
    }

    ver_r <- strsplit(as.character(R.version["version.string"]), " ")[[1]][3]
    ver_pgenlibr <- utils::packageVersion("pgenlibr")
    system(
        paste(
            "cat <<-END_VERSIONS > versions.yml",
            "\\"${task.process}\\":",
            sprintf("    r-base: %s", ver_r),
            sprintf("    pgenlibr: %s", ver_pgenlibr),
            "END_VERSIONS\\n",
            sep = "\\n"
        )
    )
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${perm_seed ? '.perm.tsv.gz' : '.tsv.gwas.gz'}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        pgenlibr: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"pgenlibr\\")))")
    END_VERSIONS
    """
}
