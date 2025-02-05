include { MANHATTAN; QQ      } from '../../modules/local/r/plots_eqtl'
include { COMBINE_EQTL_FILES } from '../../modules/local/r/combine_eqtl_files'

workflow POSTPROCESSING_EQTL {
    take:
    gwas         // channel: [mandatory] [ meta, gwas ]
    p_thr        // value  : [mandatory] [nominal p value threshold]

    main:
    versions = Channel.empty()

    // Extract only file paths and ensure they are collected as a list
    gwas.map { it[1] }  // Extract only the file path (second element)
    .collect()
    .map { file_list -> tuple(*file_list) }
    .set { gwas_file_paths }

    // Pass extracted file p aths to the COMBINE_EQTL_FILES process
    COMBINE_EQTL_FILES(gwas_file_paths)
    .set { combined_gwas }
    combined_gwas.view { "DEBUG combined_gwas: $it" }

    QQ(combined_gwas)

    MANHATTAN(combined_gwas, p_thr)

    versions.mix ( MANHATTAN.out.versions            ) .set { versions }

    emit:

    versions // channel: [ versions.yml ]
}
