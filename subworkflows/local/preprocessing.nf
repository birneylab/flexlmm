include {
    VCF_TO_PGEN; BCF_TO_PGEN; PGEN_TO_PGEN; PED_TO_PGEN; BED_TO_PGEN; BGEN_TO_PGEN
} from '../../modules/local/plink2/geno_to_pgen'
include { ESTIMATE_FREQ        } from '../../modules/local/plink2/estimate_freq'
include { GET_CHR_NAMES        } from '../../modules/local/plink2/get_chr_names'
include { SPLIT_CHR            } from '../../modules/local/plink2/split_chr'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { TRANSFORM_PHENOTYPES } from '../../modules/local/plink2/transform_phenotypes'
include { PHENO_TO_RDS         } from '../../modules/local/r/pheno_to_rds'
include { VALIDATE_FORMULAS    } from '../../modules/local/r/validate_formulas'
include { GET_DESIGN_MATRIX    } from '../../modules/local/r/get_design_matrix'
include { MATCH_SAMPLES        } from '../../modules/local/r/match_samples'


def select_chr   = params.select_chr   ? ( params.select_chr as String   ).split(",") : null
def select_pheno = params.select_pheno ? ( params.select_pheno as String ).split(",") : null
def maf_min      = params.maf_min      ?: []
def use_dosage   = params.use_dosage


workflow PREPROCESSING {
    take:
    // different genotype formats
    vcf                    // value: [ optional ]
    bcf                    // value: [ optional ]
    bgen                   // value: [ optional ]
    sample                 // value: [ optional ]
    bed                    // value: [ optional ]
    bim                    // value: [ optional ]
    fam                    // value: [ optional ]
    ped                    // value: [ optional ]
    map_f                  // value: [ optional ]
    pgen                   // value: [ optional ]
    psam                   // value: [ optional ]
    pvar                   // value: [ optional ]
    
    pheno                  // value: [mandatory] phenotypes
    covar                  // value: [optional ] covariates
    qcovar                 // value: [optional ] covariates
    freq                   // value: [optional ] vcf_file
    permute_by             // value: [optional ] permute_by factor

    null_model_formula_str // value: [mandatory] null model R formula
    model_formula_str      // value: [mandatory] model R formula

    main:
    versions = Channel.empty()
    
    // convert genotype input to pgen format
    if ( bgen && sample ){
        BGEN_TO_PGEN ( [ [id: bgen.simpleName], bgen, sample ], maf_min, use_dosage )
        BGEN_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( BGEN_TO_PGEN.out.versions ) .set { versions }
    } else if ( bcf ) {
        BCF_TO_PGEN ( [ [id: bcf.simpleName], bcf ], maf_min, use_dosage )
        BCF_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( BCF_TO_PGEN.out.versions ) .set { versions }
    } else if ( vcf ) {
        VCF_TO_PGEN ( [ [id: vcf.simpleName], vcf ], maf_min, use_dosage )
        VCF_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( VCF_TO_PGEN.out.versions ) .set { versions }
    } else if ( bed && bim && fam ) {
        BED_TO_PGEN ( [ [id: bed.simpleName], bed, bim, fam ], maf_min, use_dosage )
        BED_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( BED_TO_PGEN.out.versions ) .set { versions }
    } else if ( ped && map_f ) {
        PED_TO_PGEN ( [ [id: ped.simpleName], ped, map_f ], maf_min, use_dosage )
        PED_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( PED_TO_PGEN.out.versions ) .set { versions }
    } else if ( pgen && psam && pvar ) {
        PGEN_TO_PGEN ( [ [id: pgen.simpleName], pgen, psam, pvar ], maf_min, use_dosage )
        PGEN_TO_PGEN.out.pgen_psam_pvar.set { pgen_pvar_psam }
        versions.mix ( PGEN_TO_PGEN.out.versions ) .set { versions }
    } else {
        error (
            "No suitable combination of input genotypes has been detected. This pipeline needs one of the following combinations:\n" +
            "vcf\n" +
            "bcf\n" +
            "bgen, gen, sample\n" +
            "pgen, pvar, pvar\n" +
            "bed, bim, fam\n" +
            "ped, map (parameter is called map_f)\n"
        )
    }

    //GET_CHR_NAMES ( pgen_pvar_psam.map { it[2] } )
    //GET_CHR_NAMES.out.out
    //.map { meta, txt -> txt }
    //.splitText ()
    //.map { it.trim() }
    //.set { chr }

    //if ( select_chr ) {
    //    chr.filter { select_chr.contains ( it ) }.set { chr }
    //}

    //if ( workflow.stubRun ){
    //    Channel.of ( ["stub_chr"] ).set { chr }
    //}

    //pgen_pvar_psam.map {
    //    meta, pgen, pvar, psam ->
    //    def new_meta = meta.clone()
    //    new_meta.id = "full_genome"
    //    [new_meta, pgen, pvar, psam, []] }
    //.set { full_genome_grm_in }

    //full_genome_pgen
    //.combine ( chr )
    //.map {
    //    meta, pgen, pvar, psam, chr ->
    //    def new_meta = meta.clone()
    //    new_meta.chr = chr
    //    new_meta.id = "chr_${chr}"
    //    [new_meta, pgen, pvar, psam, chr]
    //}
    //.set { loco_grm_in }

    //if ( freq ) {
    //    ESTIMATE_FREQ ( [ [id: freq.simpleName], freq ] )
    //    ESTIMATE_FREQ.out.freq
    //    .set { freq }
    //    FULL_GRM ( full_genome_grm_in, freq )
    //    LOCO_GRM ( loco_grm_in       , freq )

    //    versions.mix ( ESTIMATE_FREQ.out.versions        ) .set { versions }
    //} else {
    //    FULL_GRM ( full_genome_grm_in, [[id:null], []] )
    //    LOCO_GRM ( loco_grm_in       , [[id:null], []] )
    //}

    //TRANSFORM_PHENOTYPES ( full_genome_pgen.combine ( [ pheno ] ) )
    //TRANSFORM_PHENOTYPES.out.pheno
    //.ifEmpty ( [ [id: pheno.simpleName], pheno ] )
    //.set { pheno }

    //PHENO_TO_RDS ( pheno )
    //PHENO_TO_RDS.out.pheno_names
    //.map { meta, pheno_names -> pheno_names }
    //.splitCsv ( header: false )
    //.flatten ()
    //.set { pheno_names }

    //if ( select_pheno ) {
    //    pheno_names.filter { select_pheno.contains ( it ) }.set { pheno_names }
    //}

    //VALIDATE_FORMULAS ( null_model_formula_str, model_formula_str )
    //VALIDATE_FORMULAS.out.fixed_effects.set { fixed_effects_formula }
    //VALIDATE_FORMULAS.out.covariates   .set { covariate_formula     }
    //VALIDATE_FORMULAS.out.intercepts   .set { intercepts            }

    //PHENO_TO_RDS.out.pheno.combine ( pheno_names ).set { pheno }
    //GET_DESIGN_MATRIX (
    //    [
    //    [id: "covar_${covar.simpleName ?: 'null'}_qcovar_${qcovar.simpleName ?: 'null'}"],
    //    covar,
    //    qcovar
    //    ],
    //    PHENO_TO_RDS.out.pheno.map { meta, pheno -> pheno }.first(),
    //    covariate_formula,
    //    fixed_effects_formula,
    //    permute_by
    //)

    //LOCO_GRM.out.grm
    //.combine ( pheno )
    //.map {
    //    meta, grm_bin, grm_id, meta2, pheno, pheno_name ->
    //    new_meta = meta.clone()
    //    new_meta.pheno = pheno_name
    //    new_meta.id = "${meta.id}_${pheno_name}"
    //    [new_meta, grm_bin, grm_id, pheno, pheno_name]
    //}
    //.combine ( GET_DESIGN_MATRIX.out.C         .map { meta, C    -> C    } )
    //.combine ( GET_DESIGN_MATRIX.out.gxe_frame .map { meta, gxe  -> gxe  } )
    //.combine ( GET_DESIGN_MATRIX.out.perm_group.map { meta, perm -> perm } )
    //.set { match_samples_in }
    //MATCH_SAMPLES ( match_samples_in )
    //MATCH_SAMPLES.out.model_terms.set { model_terms }
    //MATCH_SAMPLES.out.gxe_frame  .set { gxe_frame   }
    //MATCH_SAMPLES.out.perm_group .set { perm_group  }

    //full_genome_pgen.combine ( MATCH_SAMPLES.out.sample_ids )
    //.map {
    //    meta, pgen, psam, pvar, meta2, sample_ids ->
    //    [ meta2, pgen, psam, pvar, sample_ids, meta2.chr ]
    //}
    //.set { split_chr_in }
    //SPLIT_CHR ( split_chr_in )
    //SPLIT_CHR.out.pgen.set { chr_pheno_pgen }

    //FULL_GRM.out.grm.mix ( LOCO_GRM.out.grm ).set { all_grms }

    // Gather versions of all tools used
    //versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    //versions.mix ( FULL_GRM.out.versions             ) .set { versions }
    //versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    //versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }
    //versions.mix ( PHENO_TO_RDS.out.versions         ) .set { versions }
    //versions.mix ( VALIDATE_FORMULAS.out.versions    ) .set { versions }
    //versions.mix ( GET_DESIGN_MATRIX.out.versions    ) .set { versions }
    //versions.mix ( MATCH_SAMPLES.out.versions        ) .set { versions }
    //versions.mix ( SPLIT_CHR.out.versions            ) .set { versions }

    emit:
    //chr_pheno_pgen        // channel: [ meta, pgen, psam, pvar ]
    //model_terms           // channel: [ meta, K, y, C ]
    //gxe_frame             // channel: [ meta, gxe_frame ]
    //perm_group            // channel: [ meta, perm_group ]
    //all_grms              // channel: [ meta, grm, grm_id ]

    //fixed_effects_formula // channel: formula_rds
    //covariate_formula     // channel: formula_rds
    //intercepts            // channel: intercepts_rds


    versions              // channel: [ versions.yml ]
}
