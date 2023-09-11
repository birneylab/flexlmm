include { AIREML               } from '../../modules/local/r/aireml'
include { CHOLESKY             } from '../../modules/local/r/cholesky'
include { DECORRELATE          } from '../../modules/local/r/decorrelate'
include { FIT_MODEL            } from '../../modules/local/r/fit_model'


workflow LMM {
    take:
    chr_pgen           // channel: [mandatory] [ meta, pgen, psam, pvar ]
    model_terms        // channel: [mandatory] [ meta, K, y, C ]
    null_model_formula // value:   [mandatory] null model R formula
    model_formula      // value:   [mandatory] model R formula

    main:
    versions = Channel.empty()

    AIREML ( model_terms )

    model_terms
    .join ( AIREML.out.hsq, failOnMismatch: true, failOnDuplicate: true )
    .set { cholesky_in }
    CHOLESKY ( cholesky_in )

    model_terms
    .map { meta, K, y, C -> [meta, y, C] }
    .join ( CHOLESKY.out.chol_L, failOnMismatch: true, failOnDuplicate: true )
    .set { decorrelate_in }
    DECORRELATE ( decorrelate_in )

    decorrelate_in
    .map{ meta, y, C, L -> [meta.chr, meta, y, C, L] }
    .join(
        chr_pgen.map{ meta, pgen, psam, pvar -> [meta.chr, meta, pgen, psam, pvar] },
        failOnMismatch: true,
        failOnDuplicate: true
    )
    .map {
        chr, meta, y, C, L, meta2, pgen, psam, pvar -> [meta, y, C, L, pgen, psam, pvar]
    }
    .set { fit_model_in }
    FIT_MODEL ( fit_model_in, null_model_formula, model_formula )

    // Gather versions of all tools used
    versions.mix ( AIREML.out.versions      ) .set { versions }
    versions.mix ( CHOLESKY.out.versions    ) .set { versions }
    versions.mix ( DECORRELATE.out.versions ) .set { versions }
    versions.mix ( FIT_MODEL.out.versions   ) .set { versions }

    emit:
    //gwas = FIT_MODEL.out.gwas // channel: [ meta, gwas ]

    versions                  // channel: [ versions.yml ]
}
