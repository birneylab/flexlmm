//
// Run imputation using STITCH in a two-step process that caches imput generation
//

include { STITCH as STITCH_GENERATEINPUTS         } from '../../modules/nf-core/stitch'
include { STITCH as STITCH_IMPUTATION             } from '../../modules/nf-core/stitch'

workflow STITCH_TWO_STEP_IMPUTE {
    take:
    stitch_input      // channel: [mandatory] [ meta, posfile, chromosome_name, K, nGen ]
    collected_samples // channel: [mandatory] [ meta, collected_crams, collected_crais, cramlist ]
    reference         // channel: [mandatory] [ meta, fasta, fasta_fai ]
    random_seed       // integer: [optional ] random_seed

    main:
    versions = Channel.empty()

    def random_seed = random_seed ?: []

    stitch_input
    .map{
        meta, posfile, chromosome_name, K, nGen ->
        [meta, posfile, [], [], chromosome_name, 1, 1]
    }
    .unique()
    .set { stitch_input_1 }

    STITCH_GENERATEINPUTS ( stitch_input_1, collected_samples, reference, random_seed )

    stitch_input
    .join ( STITCH_GENERATEINPUTS.out.input, failOnMismatch: true, failOnDuplicate: true )
    .join ( STITCH_GENERATEINPUTS.out.rdata, failOnMismatch: true, failOnDuplicate: true )
    .map {
        meta, posfile, chromosome_name, K, nGen, input, rdata ->
        [
            meta,
            posfile,
            input,
            rdata,
            chromosome_name,
            K,
            nGen,
        ]
    }
    .set { stitch_input_2 }

    STITCH_IMPUTATION( stitch_input_2, [null, [], [], []], [null, [], []], random_seed )

    versions.mix ( STITCH_GENERATEINPUTS.out.versions ).set { versions }
    versions.mix ( STITCH_IMPUTATION.out.versions     ).set { versions }

    emit:
    input = STITCH_IMPUTATION.out.input
    rdata = STITCH_IMPUTATION.out.rdata
    plots = STITCH_IMPUTATION.out.plots
    vcf   = STITCH_IMPUTATION.out.vcf
    bgen  = STITCH_IMPUTATION.out.bgen

    versions // channel: [ versions.yml ]
}
