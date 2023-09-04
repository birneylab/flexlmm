//
// Run a parameter grid search
//

include { SPLIT_POSFILE                           } from '../../subworkflows/local/split_stitch_posfile'
include { STITCH_TWO_STEP_IMPUTE                  } from '../../subworkflows/local/stitch_two_step_impute'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_STITCH } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_JOINT  } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_CONCAT                         } from '../../modules/nf-core/bcftools/concat'

def random_seed = params.random_seed ?: []

workflow GRID_SEARCH {
    take:
    collected_samples // channel: [mandatory] [ meta, collected_crams, collected_crais, cramlist ]
    reference         // channel: [mandatory] [ meta, fasta, fasta_fai ]
    stitch_posfile    // channel: [mandatory] [ meta, stitch_posfile ]
    chr_list          // channel: [mandatory] list of chromosomes names

    main:
    versions = Channel.empty()

    SPLIT_POSFILE ( reference, stitch_posfile, chr_list )
    SPLIT_POSFILE.out.positions.set { positions }

    Channel.fromPath ( params.grid_search_params )
    .splitCsv( header:true )
    .set { grid_search_params }

    positions
    .combine ( grid_search_params )
    .map {
        meta, positions, chromosome_name, grid_search_params ->
        [
            [
                "id"                   : "chromosome_${chromosome_name}"                            ,
                "publish_dir_subfolder": "K_${grid_search_params.K}_nGen_${grid_search_params.nGen}",
                "params_comb"          : grid_search_params                                         ,
            ],
            positions, chromosome_name, grid_search_params.K, grid_search_params.nGen
        ]
    }
    .set { stitch_input }

    STITCH_TWO_STEP_IMPUTE ( stitch_input, collected_samples, reference, random_seed )
    STITCH_TWO_STEP_IMPUTE.out.vcf.set { stitch_vcf }
    BCFTOOLS_INDEX_STITCH ( stitch_vcf )

    stitch_vcf
    .join( BCFTOOLS_INDEX_STITCH.out.csi, failOnMismatch: true, failOnDuplicate: true )
    .map {
        meta, vcf, csi ->
        def new_meta = meta.clone()
        new_meta.id = "joint_stitch_output"
        [new_meta, vcf, csi]
    }
    .groupTuple ()
    .set { collected_vcfs }

    BCFTOOLS_CONCAT ( collected_vcfs )
    BCFTOOLS_CONCAT.out.vcf.set { genotype_vcf }
    BCFTOOLS_INDEX_JOINT( genotype_vcf )
    genotype_vcf
    .join ( BCFTOOLS_INDEX_JOINT.out.csi, failOnMismatch: true, failOnDuplicate: true )
    .set { genotype_vcf }

    versions.mix ( SPLIT_POSFILE.out.versions          ) .set { versions }
    versions.mix ( STITCH_TWO_STEP_IMPUTE.out.versions ) .set { versions }
    versions.mix ( BCFTOOLS_INDEX_STITCH.out.versions  ) .set { versions }
    versions.mix ( BCFTOOLS_CONCAT.out.versions        ) .set { versions }
    versions.mix ( BCFTOOLS_INDEX_JOINT.out.versions   ) .set { versions }

    emit:
    genotype_vcf // channel: [ meta, vcf, vcf_index ]

    versions     // channel: [ versions.yml ]

}
