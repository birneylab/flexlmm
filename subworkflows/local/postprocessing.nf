workflow POSTPROCESSING {
    take:
    gwas        // channel: [mandatory] [ meta, gwas ]
    gwas_perm   // channel: [optional ] [ meta, gwas_perm ]

    main:
    versions = Channel.empty()

    gwas.view()

    emit:

    versions // channel: [ versions.yml ]
}
