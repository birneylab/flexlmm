include { AIREML                      } from '../../modules/local/r/aireml'
include { CHOLESKY                    } from '../../modules/local/r/cholesky'
include { DECORRELATE                 } from '../../modules/local/r/decorrelate'
include { FIT_MODEL as FIT_MODEL_ORIG } from '../../modules/local/r/fit_model'
include { FIT_MODEL as FIT_MODEL_PERM } from '../../modules/local/r/fit_model'


workflow LMM {
    take:
    chr_pheno_pgen        // channel: [mandatory] [ meta, pgen, psam, pvar ]
    model_terms           // channel: [mandatory] [ meta, K, y, C ]
    gxe_frame             // channel: [mandatory] [ meta, gxe_frame ]
    perm_group            // channel: [mandatory] [ meta, perm_group ]

    fixed_effects_formula // channel: [mandatory] formula_rds

    permutation_seeds     // channel: [optional ] permutation seeds

    main:
    versions = Channel.empty()

    AIREML ( model_terms )

    model_terms
    .join ( AIREML.out.hsq,      failOnMismatch: true, failOnDuplicate: true )
    .set { cholesky_in }
    CHOLESKY ( cholesky_in )

    model_terms
    .map { meta, K, y, C -> [meta, y, C] }
    .join ( CHOLESKY.out.chol_L, failOnMismatch: true, failOnDuplicate: true )
    .set { decorrelate_in }
    DECORRELATE ( decorrelate_in )

    DECORRELATE.out.mm_rotation
    .join ( CHOLESKY.out.chol_L, failOnMismatch: true, failOnDuplicate: true )
    .join ( gxe_frame,           failOnMismatch: true, failOnDuplicate: true )
    .join ( perm_group,          failOnMismatch: true, failOnDuplicate: true )
    .join ( chr_pheno_pgen,      failOnMismatch: true, failOnDuplicate: true )
    .map {
        meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar ->
        [meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar, []]
    }
    .set { fit_model_in }
    FIT_MODEL_ORIG ( fit_model_in, fixed_effects_formula )
    FIT_MODEL_ORIG.out.gwas.set { gwas }

    fit_model_in
    .combine ( permutation_seeds )
    .map {
        meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar, fake_seed, seed ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_perm${seed}"
        new_meta.seed = seed
        new_meta.is_perm = true
        [new_meta, y, C, L, gxe_frame, perm_group, pgen, psam, pvar, seed]
    }
    .set { fit_model_perm_in }
    FIT_MODEL_PERM (fit_model_perm_in, fixed_effects_formula)
    FIT_MODEL_PERM.out.gwas.set { gwas_perm }

    // Gather versions of all tools used
    versions.mix ( AIREML.out.versions         ) .set { versions }
    versions.mix ( CHOLESKY.out.versions       ) .set { versions }
    versions.mix ( DECORRELATE.out.versions    ) .set { versions }
    versions.mix ( FIT_MODEL_ORIG.out.versions ) .set { versions }
    versions.mix ( FIT_MODEL_PERM.out.versions ) .set { versions }

    emit:
    gwas      // channel: [ meta, gwas ]
    gwas_perm // channel: [ meta, gwas_perm ]

    versions // channel: [ versions.yml ]
}
