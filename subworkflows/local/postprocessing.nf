include { GET_MIN_P_DISTRIBUTION } from '../../modules/local/r/get_min_p_distribution'


workflow POSTPROCESSING {
    take:
    gwas        // channel: [mandatory] [ meta, gwas ]
    gwas_perm   // channel: [optional ] [ meta, gwas_perm ]

    nperms      // value  : [mandatory] [number of permutations]

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

    emit:

    versions // channel: [ versions.yml ]
}
