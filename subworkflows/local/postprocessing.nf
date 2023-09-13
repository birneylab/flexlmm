include { GET_MIN_P_DISTRIBUTION } from '../../modules/local/r/get_min_p_distribution'
include { MANHATTAN; QQ          } from '../../modules/local/r/make_plots'


workflow POSTPROCESSING {
    take:
    gwas        // channel: [mandatory] [ meta, gwas ]
    gwas_perm   // channel: [mandatory] [ meta, gwas_perm ]

    nperms      // value  : [mandatory] [number of permutations]
    p_thr       // value  : [mandatory] [nominal p value threshold]

    main:
    versions = Channel.empty()

    gwas_perm
    .map {
        meta, gwas ->
        new_meta = [id: meta.pheno, pheno: meta.pheno]
        [new_meta, gwas]
    }
    .groupTuple ()
    .set { grouped_perms }
    GET_MIN_P_DISTRIBUTION ( grouped_perms, nperms )

    gwas
    .map {
        meta, gwas ->
        new_meta = [id: meta.pheno, pheno: meta.pheno]
        [new_meta, gwas]
    }
    .groupTuple ()
    .set { grouped_gwas }
    QQ ( grouped_gwas )

    grouped_gwas
    .join ( GET_MIN_P_DISTRIBUTION.out.min_p_dist, failOnMismatch: true, failOnDuplicate: true )
    .set { manhattan_in }
    MANHATTAN ( manhattan_in, p_thr )

    versions.mix ( GET_MIN_P_DISTRIBUTION.out.versions ) .set { versions }
    versions.mix ( MANHATTAN.out.versions              ) .set { versions }

    emit:

    versions // channel: [ versions.yml ]
}
