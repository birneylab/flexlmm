# eQTL Analysis Mode

## Overview

The eQTL (expression quantitative trait loci) mode enables genome-wide association analysis between genetic variants and gene expression levels. Unlike standard GWAS which tests one phenotype across all SNPs, eQTL mode tests each gene's expression as a separate phenotype, computing SNP-to-gene associations genome-wide.

**Key differences from standard GWAS mode:**
- Processes **one gene per task** (parallelizable per gene and chromosome)
- Uses **expression phenotype files** (gene expression levels across samples)
- Automatically maps SNPs to genes within a user-defined window
- Outputs **SNP-gene associations** with effect sizes

## Quick Start

```bash
nextflow run main.nf \
  -profile birneylab \
  --eqtl \
  --pheno expression_data.tsv \
  --covar covariates.tsv \
  --pgen genotypes.pgen \
  --pvar genotypes.pvar \
  --psam genotypes.psam \
  --gtf genes.gtf \
  --window 1000000 \
  --null_model_formula "y ~ cov1 + cov2" \
  --model_formula "y ~ x + cov1 + cov2" \
  --outdir results/
```

## Input Requirements

### Expression Phenotypes (`--pheno`)

Tab-separated file with:
- First column: sample IDs (must match `psam` IID column)
- Additional columns: gene expression levels (one gene per column)
- Column names: gene identifiers (e.g., ENSG00000000001)

```
IID        ENSG00000000001  ENSG00000000002
sample_1   5.2              3.1
sample_2   4.8              3.5
```

### Genotype Annotation (`--gtf`)

GTF file for SNP-to-gene mapping. Genes are mapped to SNPs within `--window` bp of their transcription start site.

### Covariates (`--covar`, `--qcovar`)

Same format as standard GWAS mode:
- **`--covar`**: categorical covariates (factors)
- **`--qcovar`**: quantitative covariates (continuous values)

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--eqtl` | false | Enable eQTL mode |
| `--window` | 100000 | SNP-to-gene mapping window (bp from TSS) |
| `--gtf` | null | Path to GTF annotation file |
| `--select_chr` | null | Limit analysis to specific chromosome(s) (comma-separated) |
| `--select_pheno` | null | Limit analysis to specific gene(s) (comma-separated) |

See [Parameters](parameters.md) for complete reference.

## Model Specification

Use the same [formula interface](../README.md#the-formula-interface) as standard GWAS, where:
- `y` represents **gene expression levels**
- `x` represents **SNP genotypes**
- Covariates are included by name

**Example:**
```
null_model_formula: y ~ age + sex + tissue
model_formula: y ~ x + age + sex + tissue
```

This tests whether SNP genotype explains additional variance in gene expression beyond age, sex, and tissue effects.

## Output

eQTL results are in `results/lmm_eqtl/`:

**Main output file:** `*.tsv.gwas.gz` (one per gene)
- Columns: chr, pos, id, ref, alt, lrt_chisq, lrt_df, pval, beta
- `beta`: effect size (formatted as `covariate~value` pairs)

**Example row:**
```
1  1000000  rs12345  A  G  12.5  1  0.0001  x~0.32
```

**Combined results:** `combined_eqtl_results.tsv.gz`
- Merged across all genes and chromosomes
- For genome-wide eQTL plotting and interpretation

## Resources

- [Main pipeline documentation](../README.md)
- [Parameter reference](parameters.md)
- [Output documentation](output.md)
- [Formula interface guide](../README.md#the-formula-interface)
