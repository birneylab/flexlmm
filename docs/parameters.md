# birneylab/flexlmm pipeline parameters

Flexible linear mixed model framework for Genome Wide Association Studies

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `vcf` | VCF file containing the sample genotypes. Can be bgzip compressed. | `string` |  | True |  |
| `pheno` | Phenotype file in PLINK2 format <details><summary>Help</summary><small>The first columns must be either FID/IID or just IID (in which case the FID is assumed to be 0). A primary header line is required and it should begin with 'FID', '#FID', 'IID', or '#IID'). Differently from the plink2 format, additional header lines (beginning with '#', not immediately followed by 'FID'/'IID') are NOT permitted before the primary header line. <br><br>IID must match sample names in `vcf`. FID is tolerated but not used. Missing values should be specified as NA and not following plink2 conventions (i.e. -9 is NOT seen as missing).<br><br>See an [example file](../assets/test_data/tsv/pheno.tsv)</small></details>| `string` |  | True |  |
| `covar` | Categorical covariates. Same format as `pheno`. <details><summary>Help</summary><small>Columns cannot be named 'x', 'y', 'ID', '#ID', 'FID', or '#FID', or have the same name of columns in `qcovar` .<br><br>See an [example file](../assets/test_data/tsv/covar.tsv)</small></details>| `string` |  |  |  |
| `qcovar` | Quantitative covariates. Same format as `pheno`. <details><summary>Help</summary><small>Columns cannot be named 'x', 'y', 'ID', '#ID', 'FID', or '#FID', or have the same name of columns in `covar` .<br><br>See an [example file](../assets/test_data/tsv/qcovar.tsv)</small></details>| `string` |  |  |  |
| `freq` | Same format as `vcf`. Genotypes to use for estimating allele frequencies. If not provided `vcf` is used. | `string` |  |  |  |
| `select_chr` | Restrict analysis to a set of chromosomes <details><summary>Help</summary><small>Comma-separated string of chromosome names in `vcf`. If not specified all chromosomes are used. If a chromosome is specified still the rest of the genome is used to evaluate the LOCO relatedness matrix.<br><br>Example: '1,2,3'</small></details>| `string` |  |  |  |
| `select_pheno` | Restrict analysis to a set of phenotypes <details><summary>Help</summary><small>Comma-separated string of column names in `pheno` to be used as phenotypes. If not specified all the columns of `pheno` are used.<br><br>Example: 'pheno1,pheno2'</small></details>| `string` |  |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |

## Statistical parameters



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `quantile_normalise` | Should the phenotypes be quantile normalised? <details><summary>Help</summary><small>Forces phenotypes to a N(0, 1) distribution.</small></details>| `boolean` |  |  |  |
| `standardise` | Should the phenotypes be mean-centered and variance-scaled? | `boolean` |  |  |  |
| `null_model_formula` | R-style formula for the outer model that you want to use as a baseline. <details><summary>Help</summary><small>A string like 'y ~ cov1'. Here 'y' can be used to refer to the phenotype, and 'x' can be used to refer to the genotype. Column names in `covar` and `qcovar` can also be used. See https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula for more information on R formulas.</small></details>| `string` |  |  |  |
| `model_formula` | R-style formula for the nested model that includes the variable of interest <details><summary>Help</summary><small>Similar to `null_model_formula`. Must contain all the terms in `null_model_formula` plus at least an extra one. You can also include GxE terms (es. 'y ~ x + x:cov1 + cov1'), and dominance terms ('y ~ x + (x == 1) + cov1'). Arithmetic operations such as '(x == 1)' must be enclosed in parentheses.</small></details>| `string` |  |  |  |
| `permutations` | Number of permutations to be performed | `integer` | 1 |  |  |
| `permute_by` | Perform permutation within pre-defined groups <details><summary>Help</summary><small>Must be a column name in `covar`. If specified, samples are exchaged in permutations only within the levels of the specified factor.</small></details>| `string` |  |  |  |
| `p_thr` | Nominal significance threshold | `number` | 0.05 |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warinig.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details>| `boolean` |  |  | True |
