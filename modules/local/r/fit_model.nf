// checks that the models are nested and that the correct terms are used where required
process FIT_MODEL {
    tag "$meta.id"
    label 'process_low'

    container 'saulpierotti-ebi/pgenlibr@sha256:0a606298c94eae8d5f6baa76aa1234fa5e7072513615d092f169029eacee5b60'

    input:
    tuple val(meta), path(y), path(C), path(L), path(gxe_frame), path(perm_group), path(pgen), path(psam), path(pvar), val(perm_seed)
    path fixed_effects_formula
    path intercepts

    output:
    tuple val(meta), path("*.gwas.tsv.gz") , emit: gwas

    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def do_permute = perm_seed ? "TRUE" : "FALSE"
    def perm_seed  = perm_seed ?: "NULL"
    """
    #!/usr/bin/env Rscript

    y <- readRDS("${y}")
    C <- readRDS("${C}")
    L <- readRDS("${L}")
    gxe_frame <- readRDS("${gxe_frame}")
    perm_group <- readRDS("${perm_group}")

    fixed_effects_formula <- readRDS("${fixed_effects_formula}")
    intercepts <- readRDS("${intercepts}")

    if ( intercepts[["model"]] == 1 & intercepts[["null_model"]] == 1 ) {
        drop_intercept <- TRUE
    } else if ( intercepts[["model"]] == 1 & intercepts[["null_model"]] == 0 ) {
        drop_intercept <- FALSE
    } else if ( intercepts[["model"]] == 0 & intercepts[["null_model"]] == 0 ) {
        drop_intercept <- FALSE
    } else {
        stop("Intercept is in null_model but not in model, models are not nested")
    }

    psam <- read.table(
        "${psam}",
        sep = "\\t",
        header = TRUE,
        comment.char = "",
        check.names = FALSE
    )
    clean_colnames <- function(n){gsub("#", "", n)}
    colnames(psam) <- clean_colnames(colnames(psam))

    stopifnot(all(!is.null(names(y))))
    stopifnot(all(names(y) == rownames(C)))
    stopifnot(all(names(y) == rownames(L)))
    stopifnot(all(names(y) == colnames(L)))
    stopifnot(all(names(y) == rownames(gxe_frame)))
    stopifnot(all(names(y) == psam[["IID"]]))
    stopifnot(all(names(y) == names(perm_group)))
    stopifnot(
        sum(is.na(L)) + sum(is.na(C)) + sum(is.na(y)) + sum(is.na(gxe_frame)) + sum(is.na(perm_group)) == 0
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

    fit_null <- .lm.fit(x = C, y = y)
    ll_null  <- stats:::logLik.lm(fit_null)

    outname <- "${prefix}.gwas.tsv.gz"
    out_con <- gzfile(outname, "w")
    header <- "chr\\tpos\\tid\\tref\\talt\\tlrt_chisq\\tlrt_df\\tpval\\tbeta"
    writeLines(header, out_con)
    nvars <- pgenlibr::GetVariantCt(pgen)
    pb <- txtProgressBar(1, nvars, style = 3)

    t <- terms(fixed_effects_formula)
    var_promise <- attr(t, "variables")
    var_names <- attr(t, "term.labels")
    var_frame <- cbind(x = 0, gxe_frame)
    var_frame[["x"]] <- pgenlibr::Buf(pgen)

    for (i in 1:nvars) {
        setTxtProgressBar(pb, i)

        pgenlibr::ReadHardcalls(pgen, var_frame[["x"]], i)

        # df with all the math operations like (x == 1) evaluated
        curr_frame <- .External2(
            stats:::C_modelframe,
            t,
            rownames(var_frame),
            eval(var_promise, var_frame, NULL),
            var_names,
            NULL,
            NULL,
            NULL,
            na.fail
        )

        # computes contrasts and interactions
        X <- .External2(stats:::C_modelmatrix, t, curr_frame)

        # drop the intercept since it is already in C, cannot drop before model.matrix so
        # that contrast are calculated correctly
        if ( drop_intercept ) X <- subset(X, select = -`(Intercept)`)

        # names are lost in forwardsolve
        X_names <- colnames(X)

        # forwardsolve(L, X) without the wrapper
        X_mm <- .Internal(backsolve(L, X, ncol(L), FALSE, FALSE))

        if ( ${do_permute} ) X_mm <- X_mm[gt_order, , drop = FALSE]

        colnames(X_mm) <- X_names
        design_matrix <- cbind(X_mm, C)

        var_id <- pgenlibr::GetVariantId(pvar, i)
        var_info <- strsplit(var_id, '_')[[1]]
        chr <- var_info[[1]]
        pos <- var_info[[2]]
        ref <- var_info[[3]]
        alt <- var_info[[4]]

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
    touch ${prefix}.gwas.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
        pgenlibr: \$(Rscript -e "cat(as.character(utils::packageVersion(\\"pgenlibr\\")))")
    END_VERSIONS
    """
}
