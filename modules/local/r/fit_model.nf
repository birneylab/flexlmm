process FIT_MODEL {
    tag "$meta.id"
    label 'process_low'

    container 'saulpierotti-ebi/pgenlibr@sha256:0a606298c94eae8d5f6baa76aa1234fa5e7072513615d092f169029eacee5b60'

    input:
    tuple val(meta), path(mm), val(perm_seed)
    tuple val(meta2), path(pgen), path(pvar), path(psam)
    path model
    path model_frame
    path perm_group
    val use_dosage

    output:
    tuple val(meta), path("*.tsv.gwas.gz") , emit: gwas
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def do_permute = perm_seed ? "TRUE" : "FALSE"
    def perm_seed  = perm_seed ?: "NULL"
    def pgenlibr_read_func = use_dosage ? "Read" : "ReadHardcalls"
    """
    #!/usr/bin/env Rscript

    l <- readRDS("${mm}")
    y.mm <- l[["y"]]
    X.mm.null <- l[["X"]]
    L <- l[["L"]]

    perm_group <- readRDS("${perm_group}")
    model <- readRDS("${model}")
    var_idx <- readRDS("${var_idx}")

    psam <- read.table(
        "${psam}",
        sep = "\\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    clean_colnames <- function(n){gsub("#", "", n)}
    colnames(psam) <- clean_colnames(colnames(psam))

    pvar_table <- read.table(
        # header starts with # and comment line with ##,
        # by removing the first # I make the header visible
        text = sub("^#", "", readLines("${pvar}")), header = TRUE
    )[var_idx,]

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y.mm) == rownames(X.mm)))
    stopifnot(all(names(y.mm) == rownames(L)))
    stopifnot(all(names(y.mm) == colnames(L)))
    stopifnot(all(names(y.mm) == psam[["IID"]]))
    stopifnot(all(names(y.mm) == names(perm_group)))
    stopifnot(
        sum(is.na(L)) + sum(is.na(X.mm)) + sum(is.na(y.mm)) + sum(is.na(perm_group)) == 0
    )

    if ( ${do_permute} ) {
        set.seed(${perm_seed})
        gt_order <- 1:length(y)
        for ( curr_group in unique(perm_group) ) {
            gt_order[perm_group == curr_group] <- sample(
                gt_order[perm_group == curr_group], replace = FALSE
            )
        }
        set.seed(NULL)
    }

    pvar <- pgenlibr::NewPvar("${pvar}")
    pgen <- pgenlibr::NewPgen("${pgen}", pvar = pvar)

    fit_null <- .lm.fit(x = X.mm.null, y = y.mm)
    ll_null  <- stats:::logLik.lm(fit_null)
    resid_null <- resid(fit_null)

    outname <- "${prefix}.tsv.gwas.gz"
    out_con <- gzfile(outname, "w")
    header <- "chr\\tpos\\tid\\tref\\talt\\tlrt_chisq\\tlrt_df\\tpval\\tbeta"
    writeLines(header, out_con)

    nvars <- length(var_idx)
    pb <- txtProgressBar(1, nvars, style = 3)
    for (i in seq_along(var_idx)) {
        setTxtProgressBar(pb, i)
        the_var_idx <- var_idx[[i]]

        pgenlibr::${pgenlibr_read_func}(pgen, var_frame[["x"]], the_var_idx)

        # df with all the math operations like (x == 1) evaluated
        X <- model.matrix(model, data = df)

        # computes contrasts and interactions
        X <- .External2(stats:::C_modelmatrix, t, curr_frame)

        # drop the intercept since it is already in C, cannot drop before model.matrix so
        # that contrast are calculated correctly
        if ( drop_intercept ) X <- subset(X, select = -`(Intercept)`)

        # names are lost in forwardsolve
        X_names <- colnames(X)

        # forwardsolve(L, X) without the wrapper
        X_mm <- .Internal(backsolve(L, X, ncol(L), FALSE, FALSE))

        colnames(X_mm) <- X_names
        design_matrix <- cbind(X_mm, C)

        var_id <- pgenlibr::GetVariantId(pvar, the_var_id)
        chr <- pvar_table[["CHROM"]]
        pos <- pvar_table[["POS"]]
        ref <- pvar_table[["REF"]]
        alt <- pvar_table[["ALT"]]

        fit <- .lm.fit(x = design_matrix, y = y)
        beta <- paste(
            colnames(design_matrix),
            fit[["coefficients"]],
            collapse = ",",
            sep = "~"
        )
        ll_fit <- stats:::logLik.lm(fit)
        lrt_df <- attributes(ll_fit)[["df"]] - attributes(ll_null)[["df"]]
        lrt_chisq <- 2 * as.numeric(ll_fit - ll_null)
        pval <- (
            if (lrt_df > 0) pchisq(lrt_chisq, df = lrt_df, lower.tail = FALSE) else NA
        )

        lineout <- sprintf(
            "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",
            chr,
            pos,
            var_id,
            ref,
            alt,
            lrt_chisq,
            lrt_df,
            pval,
            beta
        )
        writeLines(lineout, out_con)
    }

    pgenlibr::ClosePgen(pgen)
    close(out_con)
    close(pb)

    # to make sure that the output has been written properly
    gwas <- read.table(outname, header = TRUE, sep = "\t")
    stopifnot(nrow(gwas) == nvars)
    stopifnot(ncol(gwas) == 9)

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
    touch ${prefix}.tsv.gwas.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        pgenlibr: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"pgenlibr\\")))")
    END_VERSIONS
    """
}
