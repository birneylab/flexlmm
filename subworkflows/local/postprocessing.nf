// Extract performance metrics and evaluate ground truth correlation

include { RENAME_VCF_TBI                       } from '../../modules/local/rename'
include { GLIMPSE2_CONCORDANCE as G2C_FULL     } from '../../modules/nf-core/glimpse2/concordance'
include { GLIMPSE2_CONCORDANCE as G2C_CHR      } from '../../modules/nf-core/glimpse2/concordance'
include { JOIN_RSQUARE_CHR                     } from '../../modules/local/join_rsquare_chr'
include { BCFTOOLS_QUERY as EXTRACT_INFO_SCORE } from '../../modules/nf-core/bcftools/query'
include { CAT_CAT as GZIP_INFO_SCORE           } from '../../modules/nf-core/cat/cat'

// Initialise parameters

def ac_bins       = params.glimpse2_concordance_ac_bins       ?: null
def allele_counts = params.glimpse2_concordance_allele_counts ?: null
def bins          = ac_bins || allele_counts ? null : params.glimpse2_concordance_bins
def min_val_gl    = params.glimpse2_min_val_gl
def min_val_dp    = params.glimpse2_min_val_dp

def groups        = params.glimpse2_groups  ? file ( params.glimpse2_groups  , checkIfExists: true ) : []
def samples       = params.glimpse2_samples ? file ( params.glimpse2_samples , checkIfExists: true ) : []

def glimpse2_concordance_input2 = [
    [id: null],
    groups,
    bins,
    ac_bins,
    allele_counts,
]

workflow POSTPROCESSING {
    take:
    genotype_vcf     // channel: [mandatory] [ meta, vcf, vcf_index ]
    ground_truth     // channel: [optional]  [ meta, vcf, vcf_index ]
    freq             // channel: [optional]  [ meta, vcf, vcf_index ]
    chr_list         // channel: [optional]  list of chromosomes to consider
    counter          // channel: [optional]  iteration counter (for snp_set_refinement)
    nchr             // channel: [optional]  number of chromosomes (for snp_set_refinement)


    main:
    versions = Channel.empty()
    rsquare  = Channel.empty()

    //
    // rsquare processing
    //

    if ( params.ground_truth_vcf ){

        // use the first genotype_vcf as freq_vcf if the latter is missing
        // since I use INFO/PAF it should be equivalent
        freq
        .ifEmpty( [] )
        .map { [ it ] }
        .combine ( genotype_vcf.first().map { [ it ] } )
        .map {
            freq, genotype_vcf -> freq ?: genotype_vcf
        }
        .map {
            meta, vcf, tbi ->
            def new_meta = meta.clone()
            new_meta.id = "freq"

            [ new_meta, vcf, tbi ]
        }
        .set { freq }

        // to avoid name collisions
        RENAME_VCF_TBI ( freq )
        RENAME_VCF_TBI.out.renamed.set { freq }

        genotype_vcf
        .combine ( ground_truth.map { meta, vcf, tbi -> [vcf, tbi] } )
        .combine ( freq        .map { meta, vcf, tbi -> [vcf, tbi] } )
        .map { it + [ samples ] }
        .set { glimpse2_concordance_input }

        G2C_FULL (
            glimpse2_concordance_input.combine ( chr_list.map { [ it ] } ),
            glimpse2_concordance_input2,
            min_val_gl,
            min_val_dp
        )

        G2C_FULL.out.rsquare_grp.set { rsquare_grp }
        G2C_FULL.out.rsquare_spl.set { rsquare_spl }

        // see https://github.com/odelaneau/GLIMPSE/issues/179
        glimpse2_concordance_input
        .combine ( chr_list.flatten () )
        .map {
            def chr  = it[8]
            def meta = it[0].clone()
            meta.id  = chr

            [ meta ] + it[1..8]
        }
        .set { glimpse2_concordance_input_chr }

        G2C_CHR (
            glimpse2_concordance_input_chr,
            glimpse2_concordance_input2,
            min_val_gl,
            min_val_dp
        )

        if ( params.mode != "snp_set_refinement" ) {

            G2C_CHR.out.rsquare_per_site
            .map {
                meta, rsquare_per_site ->
                def new_meta = meta.clone()
                new_meta.id = "joint_stitch_output"
                [ new_meta, rsquare_per_site ]
            }
            .groupTuple ()
            .set { join_rsquare_chr_in }

        } else {

            G2C_CHR.out.rsquare_per_site
            .map {
                meta, rsquare_per_site ->
                def new_meta = meta.clone()
                new_meta.id = "joint_stitch_output"
                [ new_meta, rsquare_per_site ]
            }
            .combine ( nchr )
            .merge   ( counter )
            // workaround since I cannot use groupTuple in recursion
            .buffer { meta, rsquare_per_site, nchr, counter -> counter == nchr }
            .map {
                buf ->
                def metas    = buf.collect { meta, rsquare, nchr, counter -> meta    }.unique ()
                def rsquares = buf.collect { meta, rsquare, nchr, counter -> rsquare }.unique ()
                def nchrs    = buf.collect { meta, rsquare, nchr, counter -> nchr    }.unique ()

                assert metas.size() == 1
                assert nchrs.size() == 1

                def meta = metas[0]
                def nchr = nchrs[0]

                assert rsquares.size() == nchr

                def ret = [meta, rsquares]

                return(ret)
            }
            .set { join_rsquare_chr_in }

        }

        JOIN_RSQUARE_CHR ( join_rsquare_chr_in )

        JOIN_RSQUARE_CHR.out.rsquare_per_site
        .join ( rsquare_spl, failOnMismatch: true, failOnDuplicate: true )
        .join ( rsquare_grp, failOnMismatch: true, failOnDuplicate: true )
        .set { rsquare }

        versions.mix ( RENAME_VCF_TBI.out.versions   ).set { versions }
        versions.mix ( G2C_FULL.out.versions         ).set { versions }
        versions.mix ( G2C_CHR.out.versions          ).set { versions }
        versions.mix ( JOIN_RSQUARE_CHR.out.versions ).set { versions }

    }

    //
    // info_score processing
    //

    Channel.of ( 'chr,pos,ref,alt,info_score' )
    .collectFile ( newLine: true )
    .set { info_score_header }

    EXTRACT_INFO_SCORE ( genotype_vcf, [], [], [] )
    EXTRACT_INFO_SCORE.out.output
    .combine ( info_score_header )
    .map {
        meta, info_score, header ->

        [ meta, [ header, info_score ] ]
    }
    .set { info_score }
    GZIP_INFO_SCORE ( info_score )
    GZIP_INFO_SCORE.out.file_out.set { info_score }

    versions.mix ( EXTRACT_INFO_SCORE.out.versions ).set { versions }
    versions.mix ( GZIP_INFO_SCORE.out.versions    ).set { versions }

    emit:
    info_score // channel: [ meta, info_score_csv ]
    rsquare    // channel: [ meta, rsquare_per_site, rsquare_spl, rsquare_grp ]

    versions   // channel: [ versions.yml ]
}
