include { COMPUTE_HERITABILITY       } from '../../modules/local/r/compute_heritability'
include { CAT_PERM; SUMMARISE_BY_CHR } from '../../modules/local/r/get_null_p_dist'
include { MANHATTAN; QQ; RELATEDNESS } from '../../modules/local/r/make_plots'


workflow POSTPROCESSING {
    take:
    gwas         // channel: [mandatory] [ meta, gwas ]
    gwas_perm    // channel: [mandatory] [ meta, gwas_perm ]
    grms         // channel: [mandatory] [ meta, grm, grm_id ]
    heritability // channel: [ meta, gaston_rds ]

    nperms       // value  : [mandatory] [number of permutations]
    p_thr        // value  : [mandatory] [nominal p value threshold]

    main:
    versions = Channel.empty()

    RELATEDNESS ( grms )
    COMPUTE_HERITABILITY ( heritability.map { meta, rds -> rds }.collect() )

    gwas_perm.map {
        meta, perm ->
        new_meta = [id: meta.pheno_name, pheno: meta.pheno_name]
        [new_meta, perm]
    }
    .groupTuple ( by : 0 )
    .set { grouped_perms }
    CAT_PERM ( grouped_perms )
    SUMMARISE_BY_CHR ( CAT_PERM.out.min_p_dist, nperms )

    gwas.map {
        meta, gwas ->
        new_meta = [id: meta.pheno_name, pheno: meta.pheno_name]
        [new_meta, gwas]
    }
    .groupTuple ()
    .set { grouped_gwas }
    QQ ( grouped_gwas )

    grouped_gwas
    .join ( SUMMARISE_BY_CHR.out.min_p_dist, failOnMismatch: true, failOnDuplicate: true )
    .set { manhattan_in }
    MANHATTAN ( manhattan_in, p_thr )

    versions.mix ( RELATEDNESS.out.versions          ) .set { versions }
    versions.mix ( COMPUTE_HERITABILITY.out.versions ) .set { versions }
    versions.mix ( CAT_PERM.out.versions             ) .set { versions }
    versions.mix ( SUMMARISE_BY_CHR.out.versions     ) .set { versions }
    versions.mix ( MANHATTAN.out.versions            ) .set { versions }

    emit:

    versions // channel: [ versions.yml ]
}
