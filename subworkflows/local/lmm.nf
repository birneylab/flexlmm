include { GET_DESIGN_MATRIX    } from '../../modules/local/r/get_design_matrix'
include { GREML                } from '../../modules/local/gcta/greml'
include { CHOLESKY             } from '../../modules/local/r/cholesky'
include { DECORRELATE_PHENO    } from '../../modules/local/r/decorrelate'
include { DECORRELATE_NULL_MAT } from '../../modules/local/r/decorrelate'


workflow LMM {
    take:
    full_genome_pgen   // channel: [mandatory] [ meta, pgen, psam, pvar ]
    loco_grm           // channel: [mandatory] [ meta, grm_bin, grm_id, grm_n ]
    null_design_matrix // channel: [mandatory] [ meta, mat ]
    pheno              // value  : [mandatory] [ meta, phenotype ]

    main:
    versions = Channel.empty()

    pheno
    .map { meta, pheno -> pheno }
    .splitCsv ( header: false, limit: 1, sep: "\t" )
    .map { it[1..(it.size()-1)] } // remove #IID col
    .first ()
    .ifEmpty ( ["stub_phenotype"] )
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

    CHOLESKY.out.chol_L
    .map {
        meta, chol ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_${meta.chr}_${meta.pheno}"
        [ new_meta, chol ]
    }
    .set { chol }

    chol.combine ( pheno.map { meta, pheno -> pheno } )
    .map {
        meta, chol, pheno ->
        def pheno_col = meta.pheno
        [ meta, chol, pheno, pheno_col ]
    }
    .set { decorrelate_pheno_in }
    DECORRELATE_PHENO ( decorrelate_pheno_in )

    chol.combine ( null_design_matrix.map { meta, mat -> mat } )
    .set { decorrelate_null_mat_in }
    DECORRELATE_NULL_MAT ( decorrelate_null_mat_in )

    full_genome_pgen.view()


    // Gather versions of all tools used
    versions.mix ( GREML.out.versions                ) .set { versions }
    versions.mix ( CHOLESKY.out.versions             ) .set { versions }
    versions.mix ( DECORRELATE_PHENO.out.versions    ) .set { versions }
    versions.mix ( DECORRELATE_NULL_MAT.out.versions ) .set { versions }

    emit:

    versions          // channel: [ versions.yml ]
}