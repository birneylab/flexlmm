include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
include { ESTIMATE_FREQ        } from '../../modules/local/plink2/estimate_freq'
include { GET_CHR_NAMES        } from '../../modules/local/plink2/get_chr_names'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { TRANSFORM_PHENOTYPES } from '../../modules/local/plink2/transform_phenotypes'
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

    GET_DESIGN_MATRIX ( [ [id: "covar"], covar, qcovar ], null_model_formula )

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions          ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    versions.mix ( ESTIMATE_FREQ.out.versions        ) .set { versions }
    versions.mix ( FULL_GRM.out.versions             ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    versions.mix ( PGEN_TO_BGEN.out.versions         ) .set { versions }
    versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }

    emit:
    full_genome_pgen                               // channel: [ meta, pgen, psam, pvar ]
    loco_grm           = LOCO_GRM.out.grm          // channel: [ meta, grm_bin, grm_id, grm_n ]
    null_design_matrix = GET_DESIGN_MATRIX.out.mat // channel: [ meta, X ]
    pheno                                          // channel: [ meta, pheno ]

    versions                                       // channel: [ versions.yml ]
}
