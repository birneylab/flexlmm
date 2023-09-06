include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
include { ESTIMATE_FREQ        } from '../../modules/local/plink2/estimate_freq'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { GET_CHR_NAMES        } from '../../modules/local/get_chr_names'


workflow PREPROCESSING {
    take:
    vcf   // value: [mandatory] vcf_file
    pheno // value: [mandatory] phenotypes
    covar // value: [optional ] covariates
    freq  // value: [optional ] vcf_file

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
    .set { chr }

    ESTIMATE_FREQ ( [ [id: "freq"], freq ] )
    ESTIMATE_FREQ.out.freq
    .ifEmpty ( [ [id: "freq"], [] ] )
    .set { freq }

    FULL_GRM ( full_genome_pgen.first(), freq, []  )
    LOCO_GRM ( full_genome_pgen.first(), freq, chr )

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions   ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions ) .set { versions }
    versions.mix ( FULL_GRM.out.versions      ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions      ) .set { versions }

    emit:

    versions          // channel: [ versions.yml ]
}
