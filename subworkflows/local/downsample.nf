// Downsample high coverage cram files to desired depth

include { SAMTOOLS_COVERAGE                             } from '../../modules/nf-core/samtools/coverage'
include { GET_DOWNSAMPLE_FACTOR                         } from '../../modules/local/getdownsamplefactor'
include { SAMTOOLS_VIEW as SAMTOOLS_DOWNSAMPLE          } from '../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index'

workflow DOWNSAMPLE {
    take:
    reads               // channel: [mandatory] [ meta, cram, crai ]
    skip_chr            // channel: [mandatory] list of chromosomes to skip
    downsample_coverage // channel: [mandatory] coverage to downsample ground truth to

    main:
    versions = Channel.empty()

    reads
    .branch {
        meta, cram, crai ->
        high_cov: meta.high_cov
        low_cov : !meta.high_cov
    }
    .set { reads }

    SAMTOOLS_COVERAGE( reads.high_cov )
    GET_DOWNSAMPLE_FACTOR( SAMTOOLS_COVERAGE.out.coverage, downsample_coverage, skip_chr )
    GET_DOWNSAMPLE_FACTOR.out.downsample_factor
    .splitText ( elem: 1 ) { meta, downsample_factor -> [meta, downsample_factor.trim()] }
    .join ( reads.high_cov, failOnMismatch: true, failOnDuplicate: true )
    .map {
        meta, downsample_factor, cram, crai ->
        new_meta = meta.clone()
        new_meta.downsample_factor = downsample_factor as Float
        new_meta.id = meta.id + ".downsample"
        [new_meta, cram, crai]
    }
    .branch {
        meta, cram, crai ->
        already_downsampled: meta.downsample_factor as Float >= 1
        others: meta.downsample_factor as Float < 1
    }
    .set { reads_to_downsample }

    SAMTOOLS_DOWNSAMPLE ( reads_to_downsample.others, [["id": null], []], [] )
    SAMTOOLS_INDEX ( SAMTOOLS_DOWNSAMPLE.out.cram )
    SAMTOOLS_DOWNSAMPLE.out.cram
    .join ( SAMTOOLS_INDEX.out.crai, failOnMismatch: true, failOnDuplicate: true )
    .set { downsampled_reads }

    reads.low_cov
    .mix ( downsampled_reads )
    .mix ( reads_to_downsample.already_downsampled )
    .set { reads }

    // Gather versions of all tools used
    versions.mix ( SAMTOOLS_COVERAGE.out.versions     ) .set { versions }
    versions.mix ( GET_DOWNSAMPLE_FACTOR.out.versions ) .set { versions }
    versions.mix ( SAMTOOLS_DOWNSAMPLE.out.versions   ) .set { versions }
    versions.mix ( SAMTOOLS_INDEX.out.versions        ) .set { versions }

    emit:
    reads    // channel: [ meta, cram, crai ]

    versions // channel: [ versions.yml ]
}
