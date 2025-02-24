include {
    VCF_TO_PGEN; BCF_TO_PGEN; PGEN_TO_PGEN; PED_TO_PGEN; BED_TO_PGEN; BGEN_TO_PGEN
} from '../../modules/local/plink2/geno_to_pgen'
include { GET_CHR_NAMES        } from '../../modules/local/plink2/get_chr_names'
include { GET_VAR_IDX_EQTL     } from '../../modules/local/r/get_var_idx_eqtl'
include { MAKE_GRM as LOCO_GRM } from '../../modules/local/plink2/make_grm'
include { MAKE_GRM as FULL_GRM } from '../../modules/local/plink2/make_grm'
include { GENE_EXPR_TO_PHENO   } from '../../modules/local/r/gene_expr_to_pheno'
include { TRANSFORM_PHENOTYPES } from '../../modules/local/plink2/transform_phenotypes'
include { PHENO_TO_RDS         } from '../../modules/local/r/pheno_to_rds'
include { VALIDATE_FORMULAS    } from '../../modules/local/r/validate_formulas'
include { GET_DESIGN_MATRIX    } from '../../modules/local/r/get_design_matrix'

def select_chr   = params.select_chr   ? ( params.select_chr as String   ).split(",") : null
def select_pheno = params.select_pheno ? ( params.select_pheno as String ).split(",") : null
def maf_min      = params.maf_min      ?: []
def use_dosage   = params.use_dosage



workflow PREPROCESSING_EQTL {
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
    freq                   // value: [optional ] plink2 freq file

    gtf               // value: [ optional ] GTF annotation file
    window            // value: [ optional ] SNP-to-gene mapping window

    null_model_formula_str // value: [mandatory] null model R formula
    model_formula_str      // value: [mandatory] model R formula

    main:
    versions = Channel.empty()
    
    // convert genotype input to pgen format
    if ( bgen && sample ){
        BGEN_TO_PGEN ( [ [id: bgen.simpleName], bgen, sample ], maf_min, use_dosage )
        BGEN_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
        versions.mix ( BGEN_TO_PGEN.out.versions ) .set { versions }
    } else if ( bcf ) {
        BCF_TO_PGEN ( [ [id: bcf.simpleName], bcf ], maf_min, use_dosage )
        BCF_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
        versions.mix ( BCF_TO_PGEN.out.versions ) .set { versions }
    } else if ( vcf ) {
        VCF_TO_PGEN ( [ [id: vcf.simpleName], vcf ], maf_min, use_dosage )
        VCF_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
        versions.mix ( VCF_TO_PGEN.out.versions ) .set { versions }
    } else if ( bed && bim && fam ) {
        BED_TO_PGEN ( [ [id: bed.simpleName], bed, bim, fam ], maf_min, use_dosage )
        BED_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
        versions.mix ( BED_TO_PGEN.out.versions ) .set { versions }
    } else if ( ped && map_f ) {
        PED_TO_PGEN ( [ [id: ped.simpleName], ped, map_f ], maf_min, use_dosage )
        PED_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
        versions.mix ( PED_TO_PGEN.out.versions ) .set { versions }
    } else if ( pgen && psam && pvar ) {
        PGEN_TO_PGEN ( [ [id: pgen.simpleName], pgen, psam, pvar ], maf_min, use_dosage )
        PGEN_TO_PGEN.out.pgen_pvar_psam.set { pgen_pvar_psam }
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

    // get the names of the chromosomes from the pvar file
    GET_CHR_NAMES ( pgen_pvar_psam.map { meta, pgen, pvar, psam -> [ meta, pvar ] } )
    GET_CHR_NAMES.out.out.map { meta, txt -> txt }
    .splitText ()
    .map { it.trim() }
    .set { chr }

    // optionally filter one specific chromosome
    if ( select_chr ) {
        chr.filter { select_chr.contains ( it ) }.set { chr }
    }
    
    // for a stub run create a fake chromosome
    if ( workflow.stubRun ){
        Channel.of ( ["stub_chr"] ).set { chr }
    }
    
    // change the meta id and prepare full-genome GRM input
    pgen_pvar_psam.map {
        meta, pgen, pvar, psam ->
        def new_meta = meta.clone()
        new_meta.chr = "full_genome"
        new_meta.id = "full_genome"
        [new_meta, pgen, pvar, psam, []] }
    .set { full_genome_grm_in }

    // bind the full-genome pgen with a chromosome
    // indicator and add it to the meta map for LOCO GRM calculations
    pgen_pvar_psam.combine ( chr )
    .map {
        meta, pgen, pvar, psam, chr ->
        def new_meta = meta.clone()
        new_meta.chr = chr
        new_meta.id = chr
        [new_meta, pgen, pvar, psam, chr]
    }
    .set { loco_grm_in }

    // build full genome and LOCO GRMs
    FULL_GRM ( full_genome_grm_in, freq )
    LOCO_GRM ( loco_grm_in       , freq )
    FULL_GRM.out.grm.mix ( LOCO_GRM.out.grm ).set { all_grms }
    FULL_GRM.out.grm.set { full_genome_grm }

    // convert the expression table to a PLINK2-compatible phenotype file
    GENE_EXPR_TO_PHENO( pheno )
    GENE_EXPR_TO_PHENO.out.gene_pheno.set {pheno}
    
    // standardise or quantile normalise the phenotype if requested
    TRANSFORM_PHENOTYPES ( pgen_pvar_psam.combine ( pheno ) )
    TRANSFORM_PHENOTYPES.out.pheno
    .ifEmpty ( [ [id: pheno.simpleName], pheno ] )
    .set { pheno }

    // convert phenotype to a format easy to read in R
    // and extract names of phenotypes
    PHENO_TO_RDS ( pheno )
    PHENO_TO_RDS.out.pheno_names
    .map { meta, pheno_names -> pheno_names }
    .splitCsv ( header: false )
    .flatten ()
    .set { pheno_names }

    // select a specific phenotype if requested
    if ( select_pheno ) {
        pheno_names.filter { select_pheno.contains ( it ) }.set { pheno_names }
    }
    PHENO_TO_RDS.out.pheno.combine ( pheno_names )
    .map {
        meta, pheno, pheno_name ->
        def new_meta = meta.clone()
        new_meta.pheno_name = pheno_name
        new_meta.id = pheno_name
        [ new_meta, pheno, pheno_name ]
    }
    .set { pheno }


    // check that the null model formula is nested into the model formula
    // and save them as RDS files
    VALIDATE_FORMULAS ( null_model_formula_str, model_formula_str )
    VALIDATE_FORMULAS.out.model     .set { model }
    VALIDATE_FORMULAS.out.null_model.set { null_model }

    // get the null model design matrix and the full model frame (covar + qcovar)
    GET_DESIGN_MATRIX (
        [
            [id: "covar_qcovar"],
            covar,
            qcovar
        ],
        // I just need pheno to get the sample names in case of missing covar and qcovar, so 1 is enough
        PHENO_TO_RDS.out.pheno.map { meta, pheno -> pheno }.first(),
        model,
        null_model
    )
    GET_DESIGN_MATRIX.out.x_null     .set { x_null      }
    GET_DESIGN_MATRIX.out.model_frame.set { model_frame }
     
    // get list of pgen variant indexes to test per chromosome
    pgen_pvar_psam.combine ( chr )
    .map {
        meta, pgen, pvar, psam, chr ->
        def new_meta = meta.clone()
        new_meta.id = chr
        new_meta.chr = chr
        [ new_meta, pgen, pvar, psam, chr ]
    }
    .set { get_var_idx_in }
    pheno_names_combined = pheno_names.collect()

    GET_VAR_IDX_EQTL(get_var_idx_in,pheno_names_combined,window,gtf)
    GET_VAR_IDX_EQTL.out.var_idx.set { var_idx }
   
    GET_VAR_IDX_EQTL.out.chr_pheno_map
    .filter { meta, file_path ->
        // Check if file_path is not empty beyond the header line
        def is_not_empty = file_path.size() > 0 && file_path.text.split("\n").size() > 1
        if (!is_not_empty) {
            println "Skipping chromosome ${meta.chr} due to an empty var_idx file"
        }
        return is_not_empty
    }
    .flatMap { meta, file_path ->
        def chr = meta.chr
        def genes = file_path.text.split("\n")[1..-1]
            .collect { line -> line.split("\t")[1]?.trim() }
            .findAll { it }
        genes.collect { gene -> [chr, gene] }
    }
    .set { chr_gene_pairs }

    // preprocess chr_gene_pairs into a map for quick lookup
    chr_gene_pairs
    .map { pair ->
        [pair[0], pair[1]]  
    }
    .reduce([:]) { acc, entry ->
        def (chr, gene) = entry
        acc[chr] = (acc[chr] ?: []) << gene  
        acc
    }
    .set { chr_gene_map }  // creates a map: [chr: [gene1, gene2, ...]]

    // input for aireml in lmm subworkflow
    all_grms.combine ( pheno )
    .map {
        meta1, grm_bin, grm_id, meta2, pheno, pheno_name ->
        def new_meta = meta1.clone()
        new_meta.pheno_name = pheno_name
        new_meta.id = "${meta1.chr}_${pheno_name}"
        [ new_meta, grm_bin, grm_id, pheno, pheno_name ]
    }
    .set { aireml_in }

    // filter aireml_in keeping only the right chr:gene matches
    aireml_in
    .filter { entry ->
        def meta = entry[0]
        def pheno_name = entry[4]
        def chr = meta.chr
        chr_gene_map.get().containsKey(chr) && chr_gene_map.get()[chr]?.contains(pheno_name)
    }
    .set { filtered_aireml_in }
    
    // Gather versions of all tools used
    versions.mix ( GET_CHR_NAMES.out.versions        ) .set { versions }
    versions.mix ( LOCO_GRM.out.versions             ) .set { versions }
    versions.mix ( GENE_EXPR_TO_PHENO.out.versions   ).set  { versions }
    versions.mix ( TRANSFORM_PHENOTYPES.out.versions ) .set { versions }
    versions.mix ( PHENO_TO_RDS.out.versions         ) .set { versions }
    versions.mix ( VALIDATE_FORMULAS.out.versions    ) .set { versions }
    versions.mix ( GET_DESIGN_MATRIX.out.versions    ) .set { versions }
    versions.mix ( GET_VAR_IDX_EQTL.out.versions     ) .set { versions }

    emit:
    pgen_pvar_psam        // channel: [ meta, pgen, pvar, psam ]
    x_null                // channel: [ meta, x_null ]
    model_frame           // channel: [ meta, model_frame ]
    filtered_aireml_in    // channel: [ meta, grm, grm_id ]
    model                 // channel: formula_rds
    null_model            // channel: formula_rds
    var_idx               // channel: var_idx_rds

    versions              // channel: [ versions.yml ]
}
