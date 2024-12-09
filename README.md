# ![birneylab/flexlmm](docs/images/birneylab-flexlmm_name_light.png#gh-light-mode-only) ![birneylab/flexlmm](docs/images/birneylab-flexlmm_name_dark.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**birneylab/flexlmm** is a bioinformatics pipeline that runs linear mixed models for Genome-Wide Association Studies.
Its main advantage compared to similar tools is that it is very flexible in the definition of the statistical model to be used and it performs permutations to correct for multiple testing and non-normal phenotypes.
_p_-values are evaluated with a [likelyhood ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test) between a 'null model' and a 'real model'.
Both models are specified by the user with the [R formula interface](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula).

**Disclaimer**: this pipeline uses the nf-core template but it is not part of nf-core itself.

![birneylab/flexlmm_metro_map](docs/images/birneylab_flexlmm_drawing.png)

1. Convert genotypes (from various formats, see the [parameter docs](docs/parameters.md)) to `pgen` format ([`plink2`](https://www.cog-genomics.org/plink/2.0/))
1. Compute the relatedness matrix for the whole genome and each LOCO subset ([`plink2`](https://www.cog-genomics.org/plink/2.0/))
1. Verify that the statistical model specified is nested in the null model ([`R language`](https://www.r-project.org/))
1. Estimate variance components using the null model fixed effects and the relatedness matrix ([`gaston`](https://cran.r-project.org/web/packages/gaston/index.html))
1. Compute the Cholesky decomposition of the phenotype variance-covariance matrix ([`R language`](https://www.r-project.org/))
1. Remove the covariance structure from the phenotypes and fixed effect covariates ([`R language`](https://www.r-project.org/))
1. Fit the null and complete models for each SNP, and compute a _p_-value using a likelihood ratio test ([`R language`](https://www.r-project.org/))
1. Fit the model and compute _p_-values for each permutation of the residuals ([`R language`](https://www.r-project.org/))
1. Compute the genome-wide significance threshold using the Westfall–Young minP approach ([`R language`](https://www.r-project.org/))
1. Make the final plots and outputs ([`ggplot2`](https://ggplot2.tidyverse.org/)):
   - Manhattan plot of the associations
   - Quantile-quantile plots
   - Heatmap of the relatedness matrix ([`ComplexHeatmap`](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html))
   - Full-genome heritability estimates for each phenotype
   - Empirical genome-wide null _p_-value distribution
   - GWAS association results in a tsv format that can be directly loaded in [IGV](https://igv.org/)

## Capabilities

- Permutations are reproducible since the permutation index is used as a random seed
- `Plink2` used wherever possible
- Model fitting done by loading in memory only one SNP at a time directly from the `pgen` file and using low-level internal R routines to increase performance
- Can test for dominance and interactions of the genotype with fixed covariates (for example, to test for GxE and GxG)
- Can standardize or quantile-normalize the phenotypes
- Can include quantitative and/or categorical covariates

## Input formats

- Accepted genotype input formats:
   - `vcf`
   - `bcf`
   - `bgen`
   - `pgen/psam/pvar`
   - `bed/bim/fam`
   - `ped/map`
- Phenotypes:
   - [Plink2 format](https://www.cog-genomics.org/plink/2.0/input#pheno)
- Covariates:
   - [Plink2 format](https://www.cog-genomics.org/plink/2.0/input#covar)
- Allele frequency:
   - [Plink2 format](https://www.cog-genomics.org/plink/2.0/input#read_freq)

## To be implemented

- Categorical phenotypes
- eQTLs

## The formula interface

The central concept of this pipeline is that of testing weather a certain statistical model (the 'real model') is significantly better than another statistical model (the 'null model').
The ratio of likelyhoods of the 2 models follows a $\chi^2$ distribution, from which it is possible to extract a _p_-value.
This approach is known as [Likelihood-ratio test](https://en.wikipedia.org/wiki/Likelihood-ratio_test).
A condition for this test is that the null model must be nested in the real model. That is, the real model must include all the terms present in the null model plus at least one extra term.

In this pipeline both models are supplied by the user in the form of [R formulas](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula) using the parameters `null_model_formula` and `model_formula`.
In the formulas is possible to refer to the genotype of a given SNP with the term `x`, and to the phenotype with the term `y`.
An intercept is automatically included but can be removed by adding the term `0` or `-1` to the models.
Derived quantities such as a dominance term (i.e. `x==1` for genotypes encoded as 0,1,2) can be used but must be enclosed in parenthesis.
Names of the columns in the covariate and quantitative covariate file can also be used.
Consult the [parameter documentation](docs/parameters.md) to see how these files should be formatted.

A valid real model and null model formula pair is for example:

```
model_formula: y ~ x + (x==1) + x:cov1 + cov1
null_model_formula: y ~ cov1
```

Here, `cov1` is the name of one of the columns of the file specified via `--covar` or `--qcovar`.
The _p_-value that you will obtain in this case will indicate if any of the terms `x`, `(x==1)`, and `x:cov1` have an effect significantly different from 0, when an intercept and `cov1` are accounted for.

The following formulas instead are **NOT** valid:

```
model_formula: y ~ x + (x==1) + x:cov1
null_model_formula: y ~ cov1
```

This is because the null model contains the term `cov1` which is not present in the model, and thus the formulas are not nested.

**NOTE:** the null model formula cannot contain the term `x` or functions of it. The covar and qcovar file cannot contain columns named `x` or `y`.

## Technical details

This pipeline fits a Leave-One-Chromosome-Out (LOCO) mixed model using the 'null model' terms as fixed effects and the realised relatedness matrix as correlation matrix to estimate the variance components.
The variance components are then used to estimate the variance-covariance matrix of the phenotypes.
The square root of the variance-covariance matrix is determined via [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition), and this decomposition is used to rotate the response vector and the design matrix for each 'real model' before running [Ordinary Least Squares (OLS)](https://en.wikipedia.org/wiki/Ordinary_least_squares).
Except for the fact that the fixed effect SNPs are not included in the variance components estimation (for performance reasons), fitting OLS to the rotated response and design matrix is mathematically equivalent to fitting Generalised Least Squares (i.e. fitting a mixed model) to the original response and design matrix.
This is a fairly standard approach used for example [here](https://github.com/grimmlab/permGWAS).

Permutations are run on the residuals from the null model OLS fit after the Cholesky rotation.
Given the interdependence of the realised residuals due to the loss of degrees of freedom consequent to the estimation of the fixed effect coefficients, the permutation is not actually run on the residuals themselves, but on their projection to a lower dimensional space.
A scrambled phenotype is generated from the fitted null model phenotype and the permuted residuals, and this is used to test each SNP defining the null p-value distribution.
See [Abney M (2015) Permutation Testing in the Presence of Polygenic Variation, Genetic Epidemiology, 39, 249-258](https://doi.org/10.1002/gepi.21893) and the [MVNpermute repository](https://github.com/markabney/MVNpermute) for a more detailed explanation of this approach.

Significance thresholds for a given nominal significance level are reported using [Bonferroni correction](https://en.wikipedia.org/wiki/Bonferroni_correction) and Westfall–Young permutations (see [here](https://doi.org/10.1093/bioinformatics/btac455)).
If $m$ permutations are performed, the significance threshold is set as the $t$ quantile of the empirical distribution given by the minimum p-values for each permutation (in total a set of $m$ p-values), where $t$ is the nominal significance desired.

## Integration with [birneylab/stitchimpute](https://github.com/birneylab/stitchimpute)

In order to use a genotype file obtained from the **birneylab/stitchimpute** pipeline, activate the `stitch` profile with the flag `-profile stitch`.
This correctly loads the dosage information and fills missing genotypes.

## Birneylab-specific information

For ease of use, the ideal settings for stitch for medaka samples have been specified in a profile called `medaka`.
This can be activated with the flag `-profile medaka`.
Always use this profile when working with medaka samples.
If the genotypes that you are using have been obtained with the *birneylab/stitchimpute* pipeline, you should also use the `stitch` profile.
In this case the flag to be used is `-profile medaka,stitch`

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

Since just a few files are required to run this pipeline, differently from other pipelines a samplesheet is not used.
Instead, the required files are all specified with dedicated parameters.

You can run the pipeline using:

```bash
nextflow run birneylab/flexlmm \
   -profile <docker/singularity/.../institute> \
   --vcf input.vcf.gz \
   --pheno input.pheno \
   --model_formula 'y ~ x' \
   --null_model_formula 'y ~ 1' \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

> For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](docs/parameters.md).

> **Warning**:
> It is highly recommended to use the docker or singularity profile. Some processes do not have a working conda configuration.

## Testing
To run the pipeline on minimal test data, just execute the following:

```
nextflow run birneylab/flexlmm -profile test --outdir <OUTDIR>
```
This minimal command requires Docker to be available in your system.
If you prefer, you can run the pipeline using conda or Singularity (choose one of them) instead with:

```
nextflow run birneylab/flexlmm -profile test,<conda|singularity> --outdir <OUTDIR>
```

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

> birneylab/flexlmm was originally written by Saul Pierotti.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

The main citation for `birneylab/flexlmm` is:

> Saul Pierotti, Tomas Fitzgerald, Ewan Birney
>
> **FlexLMM: a Nextflow linear mixed model framework for GWAS**
>
> Preprint on _arXiv_. 2024 Oct 02. doi: [10.48550/ARXIV.2410.01533](https://doi.org/10.48550/ARXIV.2410.01533)

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
