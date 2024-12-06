include { AIREML                      } from '../../modules/local/r/aireml'
include { CHOLESKY                    } from '../../modules/local/r/cholesky'
include { DECORRELATE                 } from '../../modules/local/r/decorrelate'
include { FIT_MODEL as FIT_MODEL_ORIG } from '../../modules/local/r/fit_model'
include { FIT_MODEL as FIT_MODEL_PERM } from '../../modules/local/r/fit_model'

def use_dosage   = params.use_dosage

workflow LMM {
    take:
    pgen_pvar_psam        // channel: [mandatory] [ meta, pgen, pvar, psam ]
    model_terms           // channel: [mandatory] [ meta, K, y, X, var_range ]
    perm_group            // channel: [mandatory] [ meta, perm_group ]
    model                 // channel: [mandatory] formula_rds
    null_model            // channel: [mandatory] formula_rds
    permutation_seeds     // channel: [optional ] permutation seeds

    main:
    versions = Channel.empty()

    // compute variance components for the random effects
    AIREML ( model_terms )

    // compute the resiudal covariance V = s2e I + s2g K from the estimated variance components
    // and decompose it with a cholesky decomposition V = L * t(L)
    // L can then be used to decorrelate the residuals
    model_terms
    .join ( AIREML.out.hsq, failOnMismatch: true, failOnDuplicate: true )
    // drop full genome GRM from downstream steps, keep only LOCO
    // I want AIREML of full genome to get overall heritability
    .filter { meta, K, y, X, hsq -> meta.chr != "full_genome" }
    .set { cholesky_in }
    CHOLESKY ( cholesky_in )

    // solve X and y for L using an efficient triangular solve to rotate the model
    // to an uncorrelated space
    model_terms
    // drop full genome GRM from downstream steps, keep only LOCO
    // I want AIREML of full genome to get overall heritability
    .filter { meta, K, y, X -> meta.chr != "full_genome" }
    .map { meta, K, y, X -> [meta, y, X] }
    .join ( CHOLESKY.out.chol_L, failOnMismatch: true, failOnDuplicate: true )
    .set { decorrelate_in }
    DECORRELATE ( decorrelate_in )

    // fit linear model to the uncorrelated data. This is equivalent to
    // fitting a mixed model to the original data.
    //DECORRELATE.out.mm_rotation
    //.join ( CHOLESKY.out.chol_L, failOnMismatch: true, failOnDuplicate: true )
    //.join ( perm_group,          failOnMismatch: true, failOnDuplicate: true )
    //.join ( pgen_pvar_psam,      failOnMismatch: true, failOnDuplicate: true )
    //.map {
    //    meta, y, X, L, perm_group, pgen, psam, pvar ->
    //    [meta, y, X, L, perm_group, pgen, psam, pvar, []]
    //}
    //.set { fit_model_in }
    //FIT_MODEL_ORIG ( fit_model_in, model, null_model, use_dosage )
    //FIT_MODEL_ORIG.out.gwas.set { gwas }

    //fit_model_in
    //.combine ( permutation_seeds )
    //.map {
    //    meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar, fake_seed, seed ->
    //    def new_meta = meta.clone()
    //    new_meta.id = "${meta.id}_perm${seed}"
    //    new_meta.seed = seed
    //    new_meta.is_perm = true
    //    [new_meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar, seed]
    //}
    //.set { fit_model_perm_in }
    //FIT_MODEL_PERM ( fit_model_perm_in, fixed_effects_formula, intercepts, use_dosage )
    //FIT_MODEL_PERM.out.gwas.set { gwas_perm }

    //// Gather versions of all tools used
    //versions.mix ( AIREML.out.versions         ) .set { versions }
    //versions.mix ( CHOLESKY.out.versions       ) .set { versions }
    //versions.mix ( DECORRELATE.out.versions    ) .set { versions }
    //versions.mix ( FIT_MODEL_ORIG.out.versions ) .set { versions }
    //versions.mix ( FIT_MODEL_PERM.out.versions ) .set { versions }

    //emit:
    //gwas      // channel: [ meta, gwas ]
    //gwas_perm // channel: [ meta, gwas_perm ]

    //versions // channel: [ versions.yml ]
}
