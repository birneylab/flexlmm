# birneylab/flexlmm: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [December 10th 2024]

Major change and simplification in workflow logic with the aim of facilitating future integration of an eQTL branch.
Change of the permutation approach to be performed on the residuals instead than on the genotypes.
Deprecation of the `permute_by` parameter as grouped permutation are not possible on the low-dimensional projection of the residuals.
