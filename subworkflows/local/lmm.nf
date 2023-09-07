include { GET_DESIGN_MATRIX } from '../../modules/local/r/get_design_matrix'
include { GREML             } from '../../modules/local/gcta/greml'
include { CHOLESKY          } from '../../modules/local/r/cholesky'


workflow LMM {
    take:
    loco_grm           // channel: [mandatory] [ meta, grm_bin, grm_id ]
    null_design_matrix // channel: [mandatory] [ meta, X ]
    pheno              // value  : [mandatory] [ meta, phenotype ]

    main:
    versions = Channel.empty()

    pheno
    .map { meta, pheno -> pheno }
    .splitCsv ( header: false, limit: 1, sep: "\t" )
    .map { it[1..(it.size()-1)] } // remove #IID col
    .first ()
    .set { pheno_names }

    pheno_names
    .map { it.size() }
    .flatMap { n_phenos -> 0..(n_phenos-1) }
    .set { pheno_idx }

    loco_grm.combine ( pheno_names.map { [ it ] } )
    .combine ( pheno_idx )
    .map {
        meta, grm, grm_id, grm_n, pheno_names, pheno_idx ->
        def new_meta = meta.clone()
        new_meta.pheno = pheno_names[pheno_idx]
        [new_meta, grm, grm_id, grm_n, pheno_idx]
    }
    .set { greml_in }

    GREML (
        greml_in,
        pheno.first(),
        null_design_matrix
    )

    greml_in.map { meta, grm, grm_id, grm_n, pheno_idx -> [meta, grm, grm_id, grm_n] }
    .join ( GREML.out.hsq, failOnMismatch: true, failOnDuplicate: true )
    .set { variance_components }

    CHOLESKY ( variance_components )

    // Gather versions of all tools used
    versions.mix ( GREML.out.versions ) .set { versions }

    emit:

    versions          // channel: [ versions.yml ]
}
