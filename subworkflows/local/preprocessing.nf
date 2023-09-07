include { VCF_TO_PGEN          } from '../../modules/local/plink2/vcf_to_pgen'
include { ESTIMATE_FREQ        } from '../../modules/local/plink2/estimate_freq'
include { GET_CHR_NAMES        } from '../../modules/local/plink2/get_chr_names'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { TRANSFORM_PHENOTYPES } from '../../modules/local/plink2/transform_phenotypes'
include { GET_DESIGN_MATRIX    } from '../../modules/local/r/get_design_matrix'
include { GREML                } from '../../modules/local/gcta/greml'


workflow PREPROCESSING {
    take:
    vcf                // value: [mandatory] vcf_file
    pheno              // value: [mandatory] phenotypes
    null_model_formula // value: [mandatory] null model R formula
    covar              // value: [optional ] covariates
    qcovar             // value: [optional ] covariates
    freq               // value: [optional ] vcf_file

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
    .ifEmpty ( [ [id: "input"], pheno ] )
    .set { pheno }

    pheno
    .map { meta, pheno -> pheno }
    .splitCsv ( header: false, limit: 1, sep: "\t" )
    .map { it[1..(it.size()-1)] } // remove #IID col
    .first ()
    .set { pheno_names }

    pheno_names
    .map { it.size() }
    .flatMap { n_phenos -> 0..(n_phenos-1) }
    .set { pheno_idx }

    LOCO_GRM.out.grm
    .combine ( pheno_names.map { [ it ] } )
    .combine ( pheno_idx )
    .map {
        meta, grm, grm_id, pheno_names, pheno_idx ->
        def new_meta = meta.clone()
        new_meta.pheno = pheno_names[pheno_idx]
        [new_meta, grm, grm_id, pheno_idx]
    }
    .set { greml_in }

    GET_DESIGN_MATRIX ( [ [id: "covar"], covar, qcovar ], null_model_formula )

    GREML (
        greml_in,
        pheno,
        GET_DESIGN_MATRIX.out.gcta_qcovar
    )

    // Gather versions of all tools used
    versions.mix ( VCF_TO_PGEN.out.versions          ) .set { versions }
    versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    versions.mix ( ESTIMATE_FREQ.out.versions        ) .set { versions }
    versions.mix ( FULL_GRM.out.versions             ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }
    versions.mix ( GET_DESIGN_MATRIX.out.versions    ) .set { versions }

    emit:
    //chr_pgen
    pheno
    covar

    versions          // channel: [ versions.yml ]
}
