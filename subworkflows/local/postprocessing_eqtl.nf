include { COMBINE_EQTL_FILES } from '../../modules/local/r/combine_eqtl_files'
include { MANHATTAN; QQ      } from '../../modules/local/r/plots_eqtl'

workflow POSTPROCESSING_EQTL {
    take:
    gwas         // channel: [mandatory] [ meta, gwas ]
    p_thr        // value  : [mandatory] [nominal p value threshold]

    main:
    versions = Channel.empty()

    // Extract only file paths and ensure they are collected as a list
    // and then convert list into a single space-separated string
    gwas.map { it[1] }
    .collect()
    .map { file_list -> file_list.join(" ") }
    .set { gwas_file_paths }
    
    // Ensure COMBINE_EQTL_FILES runs once per batch
    COMBINE_EQTL_FILES( gwas_file_paths )
    .set { combined_gwas }  // Final combined file

    QQ ( combined_gwas )

    MANHATTAN ( combined_gwas, p_thr )

    versions.mix ( MANHATTAN.out.versions ) .set { versions }

    emit:

    versions // channel: [ versions.yml ]
}
