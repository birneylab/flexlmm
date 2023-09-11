include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
include { ESTIMATE_FREQ        } from '../../modules/local/plink2/estimate_freq'
include { GET_CHR_NAMES        } from '../../modules/local/plink2/get_chr_names'
include { SPLIT_CHR            } from '../../modules/local/plink2/split_chr'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { TRANSFORM_PHENOTYPES } from '../../modules/local/plink2/transform_phenotypes'
include { PHENO_TO_RDS         } from '../../modules/local/r/pheno_to_rds'
include { GET_DESIGN_MATRIX    } from '../../modules/local/r/get_design_matrix'


workflow PREPROCESSING {
    take:
    vcf                // value: [mandatory] vcf_file
    pheno              // value: [mandatory] phenotypes
    covar              // value: [optional ] covariates
    qcovar             // value: [optional ] covariates
    freq               // value: [optional ] vcf_file

    null_model_formula // value: [mandatory] null model R formula

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

    full_genome_pgen.combine ( chr )
    .map {
        meta, pgen, psam, pvar, chr ->
        new_meta = meta.clone()
        new_meta.chr = chr
        [ new_meta, pgen, psam, pvar, chr ]
    }
    .set { split_chr_in }

    SPLIT_CHR ( split_chr_in )

    SPLIT_CHR.out.pgen.join(
        SPLIT_CHR.out.pvar, failOnMismatch: true, failOnDuplicate: true
    ).combine (
        full_genome_pgen.map { meta, pgen, psam, pvar -> psam }
    )
    .map { meta, pgen, pvar, psam -> [meta, pgen, psam, pvar] }
    .set { chr_pgen }

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
    .splitCsv ( header: false )
    .flatten ()
    .set { pheno_names }

    PHENO_TO_RDS.out.pheno.combine ( pheno_names ).set { pheno }
    GET_DESIGN_MATRIX (
        [ [id: "covar"], covar, qcovar],
        PHENO_TO_RDS.out.pheno.map { meta, pheno -> pheno }.first(),
        null_model_formula
    )

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions          ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    versions.mix ( SPLIT_CHR.out.versions            ) .set { versions }
    versions.mix ( ESTIMATE_FREQ.out.versions        ) .set { versions }
    versions.mix ( FULL_GRM.out.versions             ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }
    versions.mix ( PHENO_TO_RDS.out.versions         ) .set { versions }
    versions.mix ( GET_DESIGN_MATRIX.out.versions    ) .set { versions }

    emit:
    chr_pgen                                       // channel: [ meta, pgen, psam, pvar ]
    loco_grm           = LOCO_GRM.out.grm          // channel: [ meta, grm_bin, grm_id ]
    null_design_matrix = GET_DESIGN_MATRIX.out.mat // channel: [ meta, covariate_mat ]
    pheno                                          // channel: [ meta, pheno, pheno_name ]

    versions                                       // channel: [ versions.yml ]
}
