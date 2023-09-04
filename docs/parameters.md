# birneylab/stitchimpute pipeline parameters

A pipeline for imputing genotypes using STITCH, evaluating imputation performance against a ground truth, optimising imputation parameters, and refining the SNP set used

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row. See [usage docs](https://github.com/birneylab/stitchimpute/usage#samplesheet-input).</small></details>| `string` |  | True |  |
| `mode` | Which branch of the pipeline to run | `string` | imputation |  |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `skip_chr` | Chromosomes to skip <details><summary>Help</summary><small>Chromosomes from the reference provided that should not be imputed. Multiple chromosomes can be separated by commas. Es. "1,2,3".</small></details>| `string` |  |  |  |
| `grid_search_params` | CSV files containing the values of the K and nGen parameters to be used for the grid search workflow <details><summary>Help</summary><small>Must have column names `K` and `nGen`. Each line corresponds to a combination of parameters to be tested.</small></details>| `string` |  |  |  |
| `snp_filtering_criteria` | CSV file containing the threshold for inclusion of SNPs in the iterative refinement process <details><summary>Help</summary><small>TBD</small></details>| `string` |  |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Stitch options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `stitch_posfile` | Positions to run the imputation over <details><summary>Help</summary><small>Where to find file with positions to run. File is tab seperated with no header, one row per SNP, with col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 = reference base, col 4 = alternate base. Bases are capitalized. Example first row: 1<tab>1000<tab>A<tab>G<tab></small></details>| `string` |  | True |  |
| `stitch_K` | Number of ancestral haplotypes <details><summary>Help</summary><small>See STITCH documentation for more details. Required for imputation mode.</small></details>| `integer` |  |  |  |
| `stitch_nGen` | Number of generations since founding of the population <details><summary>Help</summary><small>See STITCH documentation for more details. Required for imputation mode.</small></details>| `integer` |  |  |  |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genome` | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details>| `string` |  |  |  |
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. If you don't have a BWA index available this will be generated for you automatically. Combine with `--save_reference` to save BWA index for future runs.</small></details>| `string` |  |  |  |
| `igenomes_base` | Directory / URL base for iGenomes references. | `string` | s3://ngi-igenomes/igenomes |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` |  |  | True |

## Ground truth



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `ground_truth_vcf` | VCF file with ground truth calls <details><summary>Help</summary><small>Accepted format are vcf, vcf.gz, and bcf. Sample names must be identical to the sample names in the SM tag of the cram indicated in the samplesheet. Used to calculate per-sample correlation with the imputation results. The file is re-formatted appropriately and given in input to STITCH with the --genfile flag.</small></details>| `string` |  |  |  |
| `downsample_coverage` | To what average depth should the ground truth cram files be downsampled to? <details><summary>Help</summary><small>To what average depth should the ground truth cram files be downsampled to? If not specidied no downsampling is done. Must be a numeric value.</small></details>| `number` |  |  |  |
| `correlation_imputed_dosage_type` | Type of dosage to use in calculating the correlation <details><summary>Help</summary><small>Should the correlation with the ground truth be calculated on the posterior haplotype dosage ("soft"), or on the dosage corresponding to the posterior genotype that has maximal posterior probability ("hard")? The hard dosage is an integer [0, 1, 2], the soft dosage is a real number in [0, 2].</small></details>| `string` | soft |  | True |
| `random_seed` | Random seed used for downsampling | `integer` |  |  |  |

## Other options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `filter_var` | Variable to use for filtering the SNPs <details><summary>Help</summary><small>Can be "info_score" or "pearson_r". By default, if this is not set "pearson_r" is used if ground_truth_vcf is set, "info_score" otherwise. If ground_truth_vcf is not defined this is ignored and "info_score" is used in any case.</small></details>| `string` |  |  | True |

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
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warinig.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details>| `boolean` |  |  | True |
