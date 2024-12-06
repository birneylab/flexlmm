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

WorkflowFLexlmm.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory file parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// check correct genotype input
def vcf    = params.vcf    ? file(params.vcf   , checkIfExists: true ) : null
def bcf    = params.bcf    ? file(params.bcf   , checkIfExists: true ) : null
def bgen   = params.bgen   ? file(params.bgen  , checkIfExists: true ) : null
def sample = params.sample ? file(params.sample, checkIfExists: true ) : null
def bed    = params.bed    ? file(params.bed   , checkIfExists: true ) : null
def bim    = params.bim    ? file(params.bim   , checkIfExists: true ) : null
def fam    = params.fam    ? file(params.fam   , checkIfExists: true ) : null
def ped    = params.ped    ? file(params.ped   , checkIfExists: true ) : null
def map_f  = params.map_f  ? file(params.map_f , checkIfExists: true ) : null
def pgen   = params.pgen   ? file(params.pgen  , checkIfExists: true ) : null
def psam   = params.psam   ? file(params.psam  , checkIfExists: true ) : null
def pvar   = params.pvar   ? file(params.pvar  , checkIfExists: true ) : null
if (
    !(
        (  pgen &&  psam &&  pvar && !bgen && !sample && !bed && !bim && !fam && !ped && !map_f && !vcf && !bcf) ||
        ( !pgen && !psam && !pvar &&  bgen &&  sample && !bed && !bim && !fam && !ped && !map_f && !vcf && !bcf) ||
        ( !pgen && !psam && !pvar && !bgen && !sample &&  bed &&  bim &&  fam && !ped && !map_f && !vcf && !bcf) ||
        ( !pgen && !psam && !pvar && !bgen && !sample && !bed && !bim && !fam &&  ped &&  map_f && !vcf && !bcf) ||
        ( !pgen && !psam && !pvar && !bgen && !sample && !bed && !bim && !fam && !ped && !map_f &&  vcf && !bcf) ||
        ( !pgen && !psam && !pvar && !bgen && !sample && !bed && !bim && !fam && !ped && !map_f && !vcf &&  bcf)
    )
){
    log.error (
        "No suitable combination of input genotypes has been detected. This pipeline needs one of the following combinations and no additional genotype input files specified:\n" +
        "\tvcf\n" +
        "\tbcf\n" +
        "\tbgen, gen, sample\n" +
        "\tpgen, pvar, pvar\n" +
        "\tbed, bim, fam\n" +
        "\tped, map\n"
    )
}
def pheno = file( params.pheno, checkIfExists: true )

def null_model_formula = params.null_model_formula
def model_formula      = params.model_formula

def freq_f = params.freq   ? file(params.freq  , checkIfExists: true) : []
def freq   = freq_f ? [ [ id: freq_f.simpleName ], freq_f ] : [ [ id:null ], null ]
def covar  = params.covar  ? file(params.covar , checkIfExists: true) : []
def qcovar = params.qcovar ? file(params.qcovar, checkIfExists: true) : []

def permutation_seeds = params.permutations ? 1..params.permutations : []
def nperms            = params.permutations ?: 0
def permute_by        = params.permute_by   ?: []
def p_thr             = params.p_thr

if ( params.quantile_normalise && params.standardise ) {
    error "Activating both quantile_normalise and standardise at the same time is not allowed"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PREPROCESSING  } from '../subworkflows/local/preprocessing'
include { LMM            } from '../subworkflows/local/lmm'
//include { POSTPROCESSING } from '../subworkflows/local/postprocessing'

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

workflow FLEXLMM {
    versions = Channel.empty ()

    PREPROCESSING (
        vcf,
        bcf,
        bgen,
        sample,
        bed,
        bim,
        fam,
        ped,
        map_f,
        pgen,
        psam,
        pvar,
        pheno,
        covar,
        qcovar,
        freq,
        permute_by,
        null_model_formula,
        model_formula
    )

    LMM (
        PREPROCESSING.out.pgen_pvar_psam,
        PREPROCESSING.out.model_terms,
        PREPROCESSING.out.perm_group,
        PREPROCESSING.out.model,
        PREPROCESSING.out.null_model,
        permutation_seeds
    )

    //POSTPROCESSING (
    //    LMM.out.gwas,
    //    LMM.out.gwas_perm,
    //    PREPROCESSING.out.all_grms,
    //    nperms,
    //    p_thr
    //)

    versions.mix ( PREPROCESSING.out.versions  ) .set { versions }
    //versions.mix ( LMM.out.versions            ) .set { versions }
    //versions.mix ( POSTPROCESSING.out.versions ) .set { versions }

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
