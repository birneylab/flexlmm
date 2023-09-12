process VALIDATE_FORMULAS {
    label 'process_low'

    conda "bioconda::r-base==4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    val null_model_formula
    val model_formula

    output:
    path "*.null.model.rds" , emit: null_model
    path "*.model.rds"      , emit: model
    path "*.C_model.rds"    , emit: covariate_model
    path "*.X_model.rds"    , emit: fixed_effects_model

    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "formula"
    """
    #!/usr/bin/env Rscript

    null_model <- formula("${null_model_formula}")
    model      <- formula("${model_formula}")

    null_model_lhs <- all.vars(null_model[[2]])
    null_model_rhs <- attr(terms(null_model), which = "term.labels")
    null_model_intercept <- attr(terms(null_model), which = "intercept")

    model_lhs <- all.vars(model[[2]])
    model_rhs <- attr(terms(model), which = "term.labels")
    model_intercept <- attr(terms(model), which = "intercept")

    if (
        !all(null_model_rhs %in% model_rhs) |
        !(length(model_rhs) > length(null_model_rhs)) |
        # null_model_intercept can be set only if model_intercept is also set, otherwise
        # the models are not nested.
        !(model_intercept > null_model_intercept)
    ) {
        stop("The null_model_formula must be nested in the model_formula")
    }

    if !(model_lhs == "y" & null_model_lhs == "y") {
        stop(
            "The response variable of the model_formula",
            "and null_model_formula can only be 'y'"
        )
    }

    extra_terms  <- setdiff(model_rhs, null_model_rhs)
    common_terms <- intersect(model_rhs, null_model_rhs)

    # 0 if both intercepts are set, 1 if only model_intercept is set (reverse not
    # possible)
    extra_intercept <- model_intercept - null_model_intercept
    extra_terms  <- c(extra_intercept, paste0("(", extra_terms , ")"))
    common_terms <- c(null_model_intercept, paste0("(", common_terms, ")"))

    fixed_effects <- reformulate(extra_terms)
    covariates    <- reformulate(common_terms)

    # C is used later to refer to the null design matrix, X to refer to the extra
    # fixed_effects effects matrix. Intercept not set since already in C or X if wanted
    model      <- formula("y ~ 0 + X + C")
    null_model <- formula("y ~ 0 + C")

    message("Model:", deparse(model))
    message("Null model:", deparse(null_model))
    message("Fixed effects:", deparse(fixed_effects))
    message("Covariates:", deparse(covariates))

    saveRDS(null_model, "${prefix}.null_model.rds")
    saveRDS(model, "${prefix}.model.rds")
    saveRDS(covariates, "${prefix}.C_model.rds")
    saveRDS(fixed_effects, "${prefix}.X_model.rds")

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

    stub:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "formula"
    """
    touch ${prefix}.null_model.rds
    touch ${prefix}.model.rds
    touch ${prefix}.C_model.rds
    touch ${prefix}.X_model.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(strsplit(as.character(R.version[\\"version.string\\"]), \\" \\")[[1]][3])")
    END_VERSIONS
    """
}
