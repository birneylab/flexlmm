//
// Separate STITCH posfile chromosome-wise
//

include { SEPARATEPOSITIONSBYCHR } from '../../modules/local/separatepositionsbychr'

workflow SPLIT_POSFILE {
    take:
    reference      // channel: [mandatory] [ meta, fasta, fasta_fai ]
    stitch_posfile // channel: [mandatory] [ meta, stitch_posfile ]
    chr_list       // channel: [mandatory] list of chromosomes to consider

    main:
    versions = Channel.empty()

    stitch_posfile
    .combine ( chr_list.flatten() )
    .map {
        meta, stitch_posfile, chromosome_name ->
        def new_meta = meta.clone()
        new_meta.id = "positions_${chromosome_name}"
        [new_meta, stitch_posfile, chromosome_name]
    }
    .set{ positions_to_split }

    SEPARATEPOSITIONSBYCHR( positions_to_split )

    SEPARATEPOSITIONSBYCHR.out.positions
    .join ( positions_to_split, failOnMismatch: true, failOnDuplicate: true )
    .map {
        meta, positions_chr, stitch_posfile, chromosome_name ->
        [meta, positions_chr, chromosome_name]
    }
    .set { positions }

    // Gather versions of all tools used
    versions.mix ( SEPARATEPOSITIONSBYCHR.out.versions ) .set { versions }

    emit:
    positions       // channel: [ meta, positions, chromosome_name ]

    versions        // channel: [ versions.yml ]
}
