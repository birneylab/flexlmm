include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { GET_CHR_NAMES        } from '../../modules/local/get_chr_names'


workflow PREPROCESSING {
    take:
    vcf               // value: [mandatory] vcf_file
    pheno_cov_table   // value: [mandatory] pheno_cov_table

    main:
    versions = Channel.empty()

    VCF_TO_PGEN ( [[id: "input"], vcf ] )
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

    FULL_GRM ( full_genome_pgen.first(), []  )
    LOCO_GRM ( full_genome_pgen.first(), chr )

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions   ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions ) .set { versions }
    versions.mix ( FULL_GRM.out.versions      ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions      ) .set { versions }

    emit:

    versions          // channel: [ versions.yml ]
}
