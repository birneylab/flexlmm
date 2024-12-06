include { AIREML                      } from '../../modules/local/r/aireml'
include { DECORRELATE                 } from '../../modules/local/r/decorrelate'
include { FIT_NULL_MODEL              } from '../../modules/local/r/fit_null_model'
include { FIT_MODEL as FIT_MODEL_ORIG } from '../../modules/local/r/fit_model'
include { FIT_MODEL as FIT_MODEL_PERM } from '../../modules/local/r/fit_model'

def use_dosage   = params.use_dosage

workflow LMM {
    take:
    pgen_pvar_psam        // channel: [mandatory] [ meta, pgen, pvar, psam ]
    x_null                // channel: [mandatory] [ meta, x_null ]
    aireml_in             // channel: [mandatory] [ meta, grm_bin, grm_id, pheno, pheno_name ]
    model_frame           // channel: [mandatory] [ meta, model_frame ]
    perm_group            // channel: [mandatory] [ meta, perm_group ]
    model                 // channel: [mandatory] formula_rds
    null_model            // channel: [mandatory] formula_rds
    var_idx               // channel: [mandatory] var_idx_rds
    permutation_seeds     // channel: [optional ] permutation seeds

    main:
    versions = Channel.empty()

    // compute variance components for the random effects
    AIREML ( aireml_in, x_null )

    // compute the resiudal covariance V = s2e I + s2g K from the estimated variance components
    // and decompose it with a cholesky decomposition V = L * t(L)
    // L is then used to decorrelate the residuals
    //
    // drop full genome GRM from downstream steps, keep only LOCO
    // I want AIREML of full genome to get overall heritability
    AIREML.out.aireml.filter { meta, aireml -> meta.chr != "full_genome" }.set { decorrelate_in }
    DECORRELATE ( decorrelate_in )

    // fit the null model and calculate likelihoods and de-correlated residuals
    FIT_NULL_MODEL (DECORRELATE.out.mm )

    //// fit linear model to the uncorrelated data. This is equivalent to
    //// fitting a mixed model to the original data.
    //FIT_MODEL_ORIG (
    //    DECORRELATE.out.mm.map { it + [[]] },
    //    pgen_pvar_psam,
    //    model,
    //    model_frame,
    //    perm_group,
    //    use_dosage
    //)
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
