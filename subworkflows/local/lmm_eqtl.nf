include { AIREML                      } from '../../modules/local/r/aireml'
include { DECORRELATE                 } from '../../modules/local/r/decorrelate'
include { FIT_NULL_MODEL              } from '../../modules/local/r/fit_null_model'
include { FIT_MODEL_EQTL              } from '../../modules/local/r/fit_model_eqtl'

def use_dosage   = params.use_dosage

workflow LMM_EQTL {
    take:
    pgen_pvar_psam        // channel: [mandatory] [ meta, pgen, pvar, psam ]
    x_null                // channel: [mandatory] [ meta, x_null ]
    aireml_in             // channel: [mandatory] [ meta, grm_bin, grm_id, pheno, pheno_name ]
    model_frame           // channel: [mandatory] [ meta, model_frame ]
    model                 // channel: [mandatory] formula_rds
    null_model            // channel: [mandatory] formula_rds
    var_idx               // channel: [mandatory] var_idx_rds

    main:
    versions = Channel.empty()

    // compute variance components for the random effects
    AIREML ( aireml_in, x_null )
    AIREML.out.aireml.filter { meta, aireml -> meta.chr == "full_genome" }.set { heritability }

    // compute the resiudal covariance V = s2e I + s2g K from the estimated variance components
    // and decompose it with a cholesky decomposition V = L * t(L)
    // L is then used to decorrelate the residuals
    //
    // drop full genome GRM from downstream steps, keep only LOCO
    // I want AIREML of full genome to get overall heritability
    AIREML.out.aireml.filter { meta, aireml -> meta.chr != "full_genome" }.set { decorrelate_in }
    DECORRELATE ( decorrelate_in )

    // fit the null model and calculate likelihoods, fitted values, de-correlated residuals,
    // and U1 matrix
    FIT_NULL_MODEL ( DECORRELATE.out.mm )

    // fit linear model to the uncorrelated data. This is equivalent to
    // fitting a mixed model to the original data.
    DECORRELATE.out.mm.map { meta, mm -> [ meta.chr, meta, mm ] }
    .combine ( var_idx.map { meta, rds -> [ meta.chr, rds ] }, by: 0 )
    .map { chr, meta, mm, var_idx -> [ meta, mm, var_idx ] }
    .join ( FIT_NULL_MODEL.out.null_model, failOnMismatch: true, failOnDuplicate: true )
    .set { fit_model_in }

    FIT_MODEL_EQTL (
        fit_model_in,
        pgen_pvar_psam,
        model,
        model_frame,
        use_dosage
    )
    FIT_MODEL_EQTL.out.out.set { gwas }

    // Gather versions of all tools used
    versions.mix ( AIREML.out.versions         ) .set { versions }
    versions.mix ( DECORRELATE.out.versions    ) .set { versions }
    versions.mix ( FIT_NULL_MODEL.out.versions ) .set { versions }
    versions.mix ( FIT_MODEL_EQTL.out.versions ) .set { versions }

    emit:
    gwas          // channel: [ meta, gwas ]
    heritability  // channel: [ meta, gaston_rds ]

    versions // channel: [ versions.yml ]
}
