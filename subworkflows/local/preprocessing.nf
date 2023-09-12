include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
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


workflow PREPROCESSING {
    take:
    vcf                    // value: [mandatory] vcf_file
    pheno                  // value: [mandatory] phenotypes
    covar                  // value: [optional ] covariates
    qcovar                 // value: [optional ] covariates
    freq                   // value: [optional ] vcf_file

    null_model_formula_str // value: [mandatory] null model R formula
    model_formula_str      // value: [mandatory] null model R formula

    main:
    versions = Channel.empty()

    VCF_TO_PGEN ( [ [id: "input"], vcf ] )
    VCF_TO_PGEN.out.pgen
    .join(
        VCF_TO_PGEN.out.psam, failOnMismatch: true, failOnDuplicate: true
    )
    .join(
        VCF_TO_PGEN.out.pvar, failOnMismatch: true, failOnDuplicate: true
    )
    .set { full_genome_pgen }

    GET_CHR_NAMES ( VCF_TO_PGEN.out.pvar )
    GET_CHR_NAMES.out.out
    .map { meta, txt -> txt }
    .splitText ()
    .map { it.trim() }
    .ifEmpty ( ["stub_chr1"] )
    .set { chr }

    ESTIMATE_FREQ ( [ [id: "freq"], freq ] )
    ESTIMATE_FREQ.out.freq
    .ifEmpty ( [ [id: "freq"], [] ] )
    .set { freq }

    full_genome_pgen
    .map { meta, pgen, pvar, psam -> [meta, pgen, pvar, psam, []] }
    .set { full_genome_grm_in }

    full_genome_pgen
    .combine ( chr )
    .map {
        meta, pgen, pvar, psam, chr ->
        def new_meta = meta.clone()
        new_meta.chr = chr
        [new_meta, pgen, pvar, psam, chr]
    }
    .set { loco_grm_in }

    FULL_GRM ( full_genome_grm_in, freq )
    LOCO_GRM ( loco_grm_in       , freq )

    TRANSFORM_PHENOTYPES ( full_genome_pgen.combine ( [ pheno ] ) )
    TRANSFORM_PHENOTYPES.out.pheno
    .ifEmpty ( [ [id: "pheno"], pheno ] )
    .set { pheno }

    PHENO_TO_RDS ( pheno )
    PHENO_TO_RDS.out.pheno_names
    .map { meta, pheno_names -> pheno_names }
    .splitCsv ( header: false )
    .flatten ()
    .set { pheno_names }

    VALIDATE_FORMULAS ( null_model_formula_str, model_formula_str )
    VALIDATE_FORMULAS.out.null_model   .set { null_model_formula    }
    VALIDATE_FORMULAS.out.model        .set { model_formula         }
    VALIDATE_FORMULAS.out.fixed_effects.set { fixed_effects_formula }
    VALIDATE_FORMULAS.out.covariates   .set { covariate_formula     }

    PHENO_TO_RDS.out.pheno.combine ( pheno_names ).set { pheno }
    GET_DESIGN_MATRIX (
        [ [id: "covar"], covar, qcovar],
        PHENO_TO_RDS.out.pheno.map { meta, pheno -> pheno }.first(),
        covariate_formula
    )

    LOCO_GRM.out.grm
    .combine ( pheno )
    .map {
        meta, grm_bin, grm_id, meta2, pheno, pheno_name ->
        new_meta = meta.clone()
        new_meta.pheno = pheno_name
        new_meta.id = "${meta.id}_${meta.chr}_${pheno_name}"
        [new_meta, grm_bin, grm_id, pheno, pheno_name]
    }
    .combine ( GET_DESIGN_MATRIX.out.mat.map { meta, mat -> mat } )
    .set { match_samples_in }
    MATCH_SAMPLES ( match_samples_in )
    MATCH_SAMPLES.out.model_terms.set { model_terms }

    full_genome_pgen.combine ( MATCH_SAMPLES.out.sample_ids )
    .map {
        meta, pgen, psam, pvar, meta2, sample_ids ->
        [ meta2, pgen, psam, pvar, sample_ids, meta2.chr ]
    }
    .set { split_chr_in }
    SPLIT_CHR ( split_chr_in )
    SPLIT_CHR.out.pgen.set { chr_pheno_pgen }

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions          ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    versions.mix ( ESTIMATE_FREQ.out.versions        ) .set { versions }
    versions.mix ( FULL_GRM.out.versions             ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }
    versions.mix ( PHENO_TO_RDS.out.versions         ) .set { versions }
    versions.mix ( VALIDATE_FORMULAS.out.versions    ) .set { versions }
    versions.mix ( GET_DESIGN_MATRIX.out.versions    ) .set { versions }
    versions.mix ( MATCH_SAMPLES.out.versions        ) .set { versions }
    versions.mix ( SPLIT_CHR.out.versions            ) .set { versions }

    emit:
    chr_pheno_pgen        // channel: [ meta, pgen, psam, pvar ]
    model_terms           // channel: [ meta, K, y, C ]

    null_model_formula    // channel: formula_rds
    model_formula         // channel: formula_rds
    fixed_effects_formula // channel: formula_rds
    covariate_formula     // channel: formula_rds

    versions              // channel: [ versions.yml ]
}
