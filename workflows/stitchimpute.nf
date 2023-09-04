/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Validate input parameters
WorkflowStitchimpute.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory file parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkPathParamList = [
    params.input,
    params.fasta,
    params.stitch_posfile
]

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

fasta          = params.fasta          ? Channel.fromPath(params.fasta).first()          : Channel.empty()
stitch_posfile = params.stitch_posfile ? Channel.fromPath(params.stitch_posfile).first() : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Initialise optional parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ground_truth        = params.ground_truth_vcf    ? Channel.fromPath(params.ground_truth_vcf).first() : Channel.empty()
freq_vcf            = params.freq_vcf            ? Channel.fromPath(params.freq_vcf).first()         : Channel.empty()
skip_chr            = params.skip_chr            ? params.skip_chr.split( "," )                      : []
downsample_coverage = params.downsample_coverage ?: Channel.empty()

def ground_truth_basename = params.ground_truth_vcf ? file ( params.ground_truth_vcf ).simpleName : null

if ( ground_truth_basename == "freq" ) {
    error( "To avoid name collisions, the ground_truth_vcf file cannot be called 'freq.*'." )
}

if ( ground_truth_basename == "joint_stitch_output" ) {
    error( "To avoid name collisions, the ground_truth_vcf file cannot be called 'joint_stitch_output.*'." )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check conditionally mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

switch (params.mode) {

    case "imputation":

        if (!params.stitch_K) {
            error("No value was provided for the parameter stitch_K, which is required for the imputation workflow.")
        }
        if (!params.stitch_nGen) {
            error("No value was provided for the parameter stitch_nGen, which is required for the imputation workflow.")
        }
        if (params.grid_search_params) {
            log.warn("The parameter \"grid_search_params\" is set to but will not be used in the imputation workflow. Set the \"mode\" parameter to \"grid_search\" if you want to perform a parameter search.")
        }
        if (params.snp_filtering_criteria) {
            log.warn("The parameter \"snp_filtering_criteria\" is set to but will not be used in the imputation workflow. Set the \"mode\" parameter to \"snp_set_refinement\" if you want to refine the SNP set.")
        }

        break

    case "grid_search":

        if (params.stitch_K) {
            log.warn("The parameter stitch_K is set but will not be used in the \"grid_search_params\" workflow.")
        }
        if (params.stitch_nGen) {
            log.warn("The parameter stitch_nGen is set but will not be used in the \"grid_search_params\" workflow.")
        }
        if (!params.grid_search_params) {
            error("No value was provided for the parameter \"grid_search_params\", which is required for the grid search workflow.")
        }
        if (params.snp_filtering_criteria) {
            log.warn("The parameter \"snp_filtering_criteria\" is set to but will not be used in the grid search workflow. Set the \"mode\" parameter to \"snp_set_refinement\" if you want to refine the SNP set.")
        }

        break

    case "snp_set_refinement":
        // recursion required for this workflow
        nextflow.preview.recursion = true

        if (!params.stitch_K) {
            error("No value was provided for the parameter stitch_K, which is required for the SNP set refinement workflow.")
        }
        if (!params.stitch_nGen) {
            error("No value was provided for the parameter stitch_nGen, which is required for the SNP set refinement workflow.")
        }
        if (params.grid_search_params) {
            log.warn("The parameter \"grid_search_params\" is set to but will not be used in the SNP set_refinement workflow. Set the \"mode\" parameter to \"snp_set_refinement\" if you want to refine the SNP set.")
        }
        if (!params.snp_filtering_criteria) {
            error("No value was provided for the parameter \"snp_filtering_criteria\", which is required for the snp_set_refinement workflow.")
        }

        break

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { IMPUTATION         } from '../subworkflows/local/imputation'
include { GRID_SEARCH        } from '../subworkflows/local/grid_search'
include { SNP_SET_REFINEMENT } from '../subworkflows/local/snp_set_refinement'
include { POSTPROCESSING     } from '../subworkflows/local/postprocessing'
include { PLOTTING           } from '../subworkflows/local/plotting'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow STITCHIMPUTE {
    versions = Channel.empty ()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ( file(params.input) )
    INPUT_CHECK.out.reads.set { reads }

    //
    // SUBWORKFLOW: index reference genomoe and prepare list of samples
    //
    stitch_posfile.map { [["id": null], it] }.set { stitch_posfile }

    PREPROCESSING ( reads, fasta, skip_chr, ground_truth, freq_vcf, downsample_coverage )
    PREPROCESSING.out.collected_samples.set { collected_samples }
    PREPROCESSING.out.reference        .set { reference         }
    PREPROCESSING.out.chr_list         .set { chr_list          }
    PREPROCESSING.out.ground_truth     .set { ground_truth      }
    PREPROCESSING.out.freq_vcf         .set { freq_vcf          }

    switch (params.mode) {
        case "imputation":

            //
            // SUBWORKFLOW: run the imputation
            //

            IMPUTATION ( collected_samples, reference, stitch_posfile, chr_list )
            IMPUTATION.out.genotype_vcf.set { genotype_vcf }

            ////
            //// SUBWORKFLOW: calculate ground truth correlation and make plots
            ////

            POSTPROCESSING ( genotype_vcf, ground_truth, freq_vcf, chr_list, null, null )

            POSTPROCESSING.out.rsquare   .set { rsquare    }
            POSTPROCESSING.out.info_score.set { info_score }

            versions.mix ( IMPUTATION.out.versions     ).set { versions }
            versions.mix ( POSTPROCESSING.out.versions ).set { versions }

            break

        case "grid_search":

            //
            // SUBWORKFLOW: optimise parameters
            //

            GRID_SEARCH ( collected_samples, reference, stitch_posfile, chr_list )
            GRID_SEARCH.out.genotype_vcf.set { genotype_vcf }

            //
            // SUBWORKFLOW: calculate ground truth correlation and make plots
            //

            POSTPROCESSING ( genotype_vcf, ground_truth, freq_vcf, chr_list, null, null )

            POSTPROCESSING.out.rsquare   .set { rsquare    }
            POSTPROCESSING.out.info_score.set { info_score }

            versions.mix ( GRID_SEARCH.out.versions    ).set { versions }
            versions.mix ( POSTPROCESSING.out.versions ).set { versions }

            break

        case "snp_set_refinement":

            //
            // SUBWORKFLOW: refine SNP set
            //

            SNP_SET_REFINEMENT (
                collected_samples,
                reference,
                stitch_posfile,
                chr_list,
                ground_truth,
                freq_vcf
            )
            SNP_SET_REFINEMENT.out.genotype_vcf.set { genotype_vcf }
            SNP_SET_REFINEMENT.out.performance .set { performance  }
            SNP_SET_REFINEMENT.out.info_score  .set { info_score   }
            SNP_SET_REFINEMENT.out.rsquare     .set { rsquare      }

            versions.mix ( SNP_SET_REFINEMENT.out.versions ).set { versions }

            break

    }

    PLOTTING ( info_score, rsquare )

    PLOTTING.out.plots.set { plots }

    versions.mix ( INPUT_CHECK.out.versions    ).set { versions }
    versions.mix ( PREPROCESSING.out.versions  ).set { versions }
    versions.mix ( PLOTTING.out.versions       ).set { versions }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
