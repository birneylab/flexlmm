# birneylab/flexlmm: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Convert genotypes](#convert-genotypes) - Convert genotypes to the pgen format
- [Prepare inputs](#prepare-inputs) - Format and match the data for downstream processing
- [Relatedness](#relatedness) - Computes full-genome and LOCO relatedness matrices
- [Variance components](#variance-components) - Estimate genetic and residuals variances
- [GWAS](#gwas) - Compute _p_-values and other SNP-wise statistics
- [Permutations](#permutations) - Compute the distribution of the minimum _p_-values for each permutation
- [Plots](#plots) - Plots of association results, _p_-value distributions, and relatedness matrices
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Convert genotypes

Converts the vcf genotypes in pgen format.

<details markdown="1">
<summary>Output files</summary>

- `genotypes/`
  - `*.pgen`: binary plink2 file
  - `*.psam`: sample ids
  - `*.pvar.zst`: variant information compressed with zstandard

</details>

### Prepare inputs

Matches samples between the various files, and expands categorical covariates and interaction terms. All matrices are in RDS format and can be loaded into R with the `loadRDS` function.

<details markdown="1">
<summary>Output files</summary>

- `model_matrices/{phenotype_name}/{chromosome_name}/`
  - `*.K.rds`: LOCO relatedness matrix for the current chromosome and phenotype
  - `*.C.rds`: matrix of fixed-effects for the null model
  - `*.y.rds`: phenotype vector
  - `*.gxe_frame.matched.rds`: `data.frame` object with all the covariates converted to the correct datatypes (for computing gxe terms)
  - `*.perm_group.matched.rds`: vector of groups memberships to be respected in the permutations
  - `*.sample.id`: text file with ordered sample names, one per line

</details>

### Relatedness

Computation of the genetic relatedness matrices for the full genome and with each chromosome left out (LOCO).

<details markdown="1">
<summary>Output files</summary>

- `relatedness_matrix/`
  - `loco/{left_out_chromosome_name}`: relatedness matrix evaluated on the full genome except for the named chromosome
    - `*.rel.bin`: [binary plink format](https://www.cog-genomics.org/plink/2.0/distance)
    - `*.rel.id`: sample IDs
  - `full_genome/`: relatedness matrix evaluated on the full genome
    - `*.rel.bin`: [binary plink format](https://www.cog-genomics.org/plink/2.0/distance)
    - `*.rel.id`: sample IDs

</details>

### Variance components

Estimates of the genetic and enviromental variances for each LOCO relatedness matrix and each phenotype. Can also be used to estimate heritability. RDS objects that can be loaded into R using the `loadRDS` function. They contain the whole object retured by the call to [`gaston::lmm.aireml`](https://cran.r-project.org/web/packages/gaston/index.html).

<details markdown="1">
<summary>Output files</summary>

- `variance_components/{phenotype_name}/*.hsq.rds`

</details>

### GWAS

GWAS result as compressed tab-separated file with header line, one file for each chromosome and phenotype.

<details markdown="1">
<summary>Column description</summary>

- `chr`: chromosome name
- `pos`: position of the variant
- `ref`: reference allele
- `alt`: alternate allele
- `lrt_chisq`: $\chi^2$ value of the Likelihood-ratio test
- `lrt_df`: degrees of freedom of the Likelihood-ratio test (number of extra parameters in the real model compared to the null model)
- `lrt_p`: _p_-value computed from `lrt_chisq` and `lrt_df`
- `beta`: fixed effect sizes for all the terms (including intercept and covariates)
  - This column contains several values. For each parameter, it contains the string `variable_name~beta_value`. This is separated by commas from other variable_name-beta_value pairs. Example of the content of the `beta` column for one line: `var1~0.3,var2~0.6,(Intercept)~2,x~1.2,x==1TRUE~0.9`

</details>

<details markdown="1">
<summary>Output files</summary>

- `gwas/{phenotype_name}/*.gwas.tsv.gz`

</details>

### Permutations

For each permutation, the minimum _p_-value is collected in a data frame together with the number of tests performed for that permutation and the permutation ID (which is also the random seed used). One file per phenotype will be present.

<details markdown="1">
<summary>Output files</summary>

- `permutations/*.min_p_dist.rds`: An RDS file that can be loaded in R with the `loadRDS` function. It contains a `data.frame` object with columns `permutation,min_p,n_snps`

</details>

### Plots

Manhattan plots of the GWAS results, QQ plots, and heatmaps of the relatedness matrices.

<details markdown="1">
<summary>Output files</summary>

- `plots/`
  - `manhattan/*.png`: Manhattan plots, one per phenotype
  - `qq/*.png`: quantile-quantile (qq) plots, one per phenotype
  - `relatednss_matrix`
    - `loco/{left_out_chromosome_name}/*.png`: LOCO relatedness matrix heatmaps, one per chromosome and phenotype
    - `full_genome/*.png`: full-genome relatedness matrix heatmaps, one per phenotype

</details>


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
