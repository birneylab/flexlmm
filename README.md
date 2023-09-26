# ![birneylab/flexlmm](docs/images/birneylab-flexlmm_name_light.svg#gh-light-mode-only) ![birneylab/flexlmm](docs/images/birneylab-flexlmm_name_dark.svg#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**birneylab/flexlmm** is a bioinformatics pipeline that runs linear mixed models for Genome-Wide Association Studies.
It is not particularly fast or different from other tools, but it is very flexible in the definition of the statistical model to be used and it uses permutations to correct for multiple testing and non-normal phenotypes.
P-values are evaluated with a likelyhood ratio test among a 'null model' and a 'real model'.
Both models are specified by the user with the R formula interface.

## Technical details

This pipeline fits a Leave-One-Chromosome-Out (LOCO) mixed model using the 'null model' terms as fixed effects and the realised relatedness matrix as correlation matrix to estimate the variance components.
The variance components are then used to estimate the variance-covariance matrix of the phenotypes.
The square root of the variance-covariance matrix is determined via [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition), and this decomposition is used to rotate the response vector and the design matrix for each 'real model' before running [Ordinary Least Squares (OLS)](https://en.wikipedia.org/wiki/Ordinary_least_squares).
Minus the fact that the fixed effect SNPs are not included in the variance components estimation (for performance reasons), fitting OLS to the rotated response and design matrix is mathematically equivalent to fitting Generalised Least Squares (i.e. fitting a mixed model) to the original response and design matrix.
This is a fairly standard approach used for example [here](https://github.com/grimmlab/permGWAS).

Permutations are run on the genotype vectors jointly, so that each genotype is permuted in the same way and linkage disequilibrium is maintained. Since phenotypes and covariates are not permuted, also their relationship is not altered.
Genotype permutations are performed AFTER the Cholesky rotation, so that the relatedness structure is regressed out from the correct samples and only the residuals from that are permuted.
This corresponds to a null hypothesis where the exchangeable quantities are the genotype identities when relatedness is accounted already for.
See [here](https://doi.org/10.1186/s13059-021-02354-7) for an example of this approach being used in practice.

Significance thresholds for a given nominal significance level are reported using [Bonferroni correction](https://en.wikipedia.org/wiki/Bonferroni_correction) and Westfall–Young permutations (see [here](https://doi.org/10.1093/bioinformatics/btac455)).
If $m$ permutations are permuted, the significance threshold is set as the $t$ quantile of the empirical distribution given by the minimum p-values for each permutation (in total a set of $m$ p-values), where $t$ is the nominal significance desired.

**Disclaimer**: this pipeline uses the nf-core template but it is not part of nf-core itself.

![birneylab/stitchimpute_metro_map](docs/images/birneylab-stitchimpute_metro_map.png)

<!--
**nf-core/stitchimpute** is a bioinformatics pipeline that ...
-->

1. Downsample high-coverage cram files ([`samtools`](http://www.htslib.org/doc/samtools.html); _optional_)
2. Run joint imputation with STITCH on high and low coverage cram files ([`STITCH`](https://doi.org/10.1038/ng.3594))
3. Compare imputation results to ground truth variants ([`scikit-allel`](https://scikit-allel.readthedocs.io/en/stable/) and [`anndata`](https://anndata.readthedocs.io/en/latest/); _optional_)
4. Plot the cumulative density of several per-SNP performance metrics ([`ggplot2`](https://ggplot2.tidyverse.org/)):
   - Info score
   - Pearson $r$
   - Pearson $r^2$

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,cram,crai
/path/to/sample1.cram,/path/to/sample1.cram.crai
/path/to/sample2.cram,/path/to/sample2.cram.crai
```

Each row represents a sample with its associated cram file and crai file.

Now, you can run the pipeline using:

```bash
nextflow run birneylab/stitchimpute \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

> For more details and further functionality, please refer to the [usage documentation](docs/usage.md) and the [parameter documentation](docs/parameters.md).

<!--
> TODO: add docs
> For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/stitchimpute/usage) and the [parameter documentation](https://nf-co.re/stitchimpute/parameters).
-->

## Pipeline output

<!--
To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/stitchimpute/results) tab on the nf-core website pipeline page.
-->

For more details about the output files and reports, please refer to the
[output documentation](docs/output.md).

## Credits

<!--
nf-core/stitchimpute was originally written by Saul Pierotti.
-->

> birneylab/stitchimpute was originally written by Saul Pierotti.

<!--
> We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#stitchimpute` channel](https://nfcore.slack.com/channels/stitchimpute) (you can join with [this invite](https://nf-co.re/join/slack)).
-->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/stitchimpute for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
