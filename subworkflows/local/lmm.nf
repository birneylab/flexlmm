include { AIREML               } from '../../modules/local/r/aireml'
include { CHOLESKY             } from '../../modules/local/r/cholesky'
include { DECORRELATE          } from '../../modules/local/r/decorrelate'
include { FIT_MODEL            } from '../../modules/local/r/fit_model'


workflow LMM {
    take:
    chr_pheno_pgen        // channel: [mandatory] [ meta, pgen, psam, pvar ]
    model_terms           // channel: [mandatory] [ meta, K, y, C ]

    fixed_effects_formula // channel: [mandatory] formula_rds
    model_formula         // channel: [mandatory] formula_rds
    null_model_formula         // channel: [mandatory] formula_rds

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
    .join ( chr_pheno_pgen,      failOnMismatch: true, failOnDuplicate: true )
    .set { fit_model_in }
    FIT_MODEL ( fit_model_in, fixed_effects_formula, model_formula, null_model_formula )
    FIT_MODEL.out.gwas.set { gwas }

    // Gather versions of all tools used
    versions.mix ( AIREML.out.versions      ) .set { versions }
    versions.mix ( CHOLESKY.out.versions    ) .set { versions }
    versions.mix ( DECORRELATE.out.versions ) .set { versions }
    versions.mix ( FIT_MODEL.out.versions   ) .set { versions }

    emit:
    gwas     // channel: [ meta, gwas ]

    versions // channel: [ versions.yml ]
}
