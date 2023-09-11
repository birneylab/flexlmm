process FIT_MODEL {
    tag "$meta.id"
    label 'process_low'

    container 'saulpierotti-ebi/pgenlibr@sha256:0a606298c94eae8d5f6baa76aa1234fa5e7072513615d092f169029eacee5b60'

    input:
    tuple val(meta), path(chol_L), path(pheno), path(covs), path(pgen), path(psam), path(pvar)
    val null_model_formula
    val model_formula

    output:
    tuple val(meta), path("*.gwas.tsv.gz") , emit: gwas
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    L <- readRDS("${L}")
    y <- readRDS("${y}")
    C <- readRDS("${C}")
    pvar <- pgenlibr::NewPvar("${pvar}")
    pgen <- pgenlibr::NewPgen("${pgen}", pvar = pvar)
    buf  <- pgenlibr::Buf(pgen)

    null_model <- formula(${null_model_formula})
    model <- formula(${model_formula})

    out_con <- gzfile("${prefix}.gwas.tsv.gz", "a")

    ######################################################################################
    # Sanity checks for the models
    ######################################################################################

    null_model_rhs <- attr(terms(null_model), which = "term.labels")
    null_model_lhs <- attr(terms(null_model), which = "variables")
    model_rhs <- attr(terms(model), which = "term.labels")
    model_lhs <- attr(terms(model), which = "variables")


    if (
        !all(null_model_rhs %in% model_rhs) |
        !(length(model_rhs) > length(null_model_rhs))
    ) {
        stop("The null_model_formula must be nested in the model_formula")
    }

    if (model_lhs == "y" & null_model_lhs == "y") {
        stop(
            "The response variable of the model_formula",
            "and null_model_formula can only be 'y'"
        )
    }

    extra_terms <- setdiff(model_rhs, null_model_rhs)
    if ( !all(extra_terms %in% c("x", "d")) ) {
        stop(
            "Only 'x' and 'd' terms are allowed to be present in the model_formula",
            "and absent from the null_model_formula"
        )
    }

    to_drop <- match(null_model_rhs, model_rhs)
    model <- formula(drop.terms(terms(model), to_drop))
    model <- update(model, y ~ . + C - 1) # covariates and intercept already part of C

    null_model <- formula(y ~ C)

    message("Model:", deparse(model))
    message("Null model:", deparse(null_model))

    ######################################################################################

    fit_null <- lm(null_model)
    ll_null <- logLik(fit_null)
    lrt_df <- length(extra_terms)

    process_variant <- function (i){
        setTxtProgressBar(pb, i)

        pgenlibr::ReadHardcalls(pgen, buf, i)

        tmp <- data.frame(
            x = forwardsolve(L, buf),
            d = forwardsolve(L, (buf == 1))
        )

        fit <- lm(model, data = tmp)

        ll_fit <- logLik(fit)
        lrt_chisq <- 2 * as.numeric(ll_fit - ll_null)
        p_lrt <- pchisq(lrt_chisq, df = lrt_df, lower.tail = FALSE)

        var_id <- pgenlibr::GetVariantId(pvar, i)
        var_info <- strsplit(var_id, '_')[[1]]
        chr <- var_info[[1]]
        pos <- var_info[[2]]
        ref <- var_info[[3]]
        alt <- var_info[[4]]

        lineout <- sprintf(
            "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",
            chr,
            pos,
            var_id,
            ref,
            alt,
            lrt_chisq,
            p_lrt
        )
        writeLines(lineout, out_con)

        rm(
            tmp, fit, ll_fit, lrt_chisq,
            p_lrt, var_id, var_info, chr, pos, ref, alt, lineout
        )

        return(0)
    }

    header <- "chr\\tpos\\tid\\tref\\talt\\tlrt_chisq\\tlrt_p"
    writeLines(header, out_con)
    nvars <- pgenlibr::GetVariantCt(pgen)
    pb <- txtProgressBar(1, nvars)
    ret <- sapply(1:nvars, process_variant)

    stopifnot(all(ret == 0))

    pgenlibr::ClosePgen(pgen)
    close(out_con)

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
