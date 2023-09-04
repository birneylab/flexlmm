# birneylab/stitchimpute: Usage

<!--
# nf-core/stitchimpute: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/stitchimpute/usage](https://nf-co.re/stitchimpute/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._
-->

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'

sample,cram,crai
SAMPLE1,AEG588A1_S1_L002_R1_001.cram,AEG588A1_S1_L002_R2_001.cram.crai
SAMPLE2,AEG588A1_S1_L003_R1_002.cram,AEG588A1_S1_L003_R2_002.cram.crai
SAMPLE3,AEG588A1_S1_L004_R1_003.cram,AEG588A1_S1_L004_R2_003.cram.crai
```

## Ground truth

If a ground truth is available, provide it as a VCF file with the `ground_truth_vcf` parameter. Note that sample names in the VCF must overlap with the sample names in the cram files for them to be matched correctly.

In the samplesheet, add a column `high_cov` to specify which cram files are high coverage and should be downsampled before imputation. Set the desired sequencing depth after downsampling with the `downsample_coverage` parameter.

`samplesheet.csv`:

```csv
sample,cram,crai,high_cov
/path/to/sample1.cram,/path/to/sample1.cram.crai,true
/path/to/sample2.cram,/path/to/sample2.cram.crai,false
```

And then run the pipeline with:

```bash
nextflow run birneylab/stitchimpute \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --ground_truth_vcf /path/to/ground_truth.vcf.gz \
   --downsample_coverage 0.5
```

Downsampling is done in order to have an accurate representation of the performance of the imputation for sample with a given coverage level. For this reason, it is suggested to set `downsample_coverage` to a value similar to the coverage of most samples in the analysis cohort.

## Pipeline modes

This pipeline actually consist of three related pipelines, which can be selected with the `mode` parameter.

### Imputation

The imputation workflow is selected by setting the `mode` parameter to "imputation". This is the default setting if `mode` is not specified.
The samples are imputed using STITCH with a specific set of parameters (`K` and `nGen` at least must be set, see the documentation of [STITCH](https://github.com/rwdavies/STITCH) for more details).

### Parameter optimization

The parameter optimization workflow is selected by setting the `mode` parameter to "grid_search".
In this workflow `K` and `nGen` do not need to be specified.
A set of combinations for the 2 parameters must be provided via the `grid_search_params` parameter.
This should point to a csv file containing a column with header `K` and one with header `nGen`.
Each row contains a combination of `K` and `nGen` which the pipeline will use for imputation.
The imputation will be repeated with different parameters for as many times as there are non-header lines in the file specified.

An example of `grid_search_params` file is the following:

```
K,nGen
8,1
8,2
16,2
32,2
```

### SNP set refinement

The SNP set refinement workflow is selected by setting the `mode` parameter to "snp_set_refinement".
In this mode `K` and `nGen` need to be supplied as in the "imputation" mode.
In this mode a first imputation is done, and then the imputation performance is evaluated for each SNP internally by STITCH (info_score), or externally against a ground truth.
A second imputation is then performed only using SNPs that satisfy a certain threshold for the imputation quality metric.
Then a third imputation is performed and so on.
This is the method suggested in the original STITCH paper for obtaining a reliable SNP set (https://www.nature.com/articles/ng.3594).

The parameter `snp_filtering_criteria` must be set and point to a csv file with no header, and a single column.
Each row should contain a threshold value that is used to filter out SNPs with a performance metric lower than the threshold.
The metric considered is the info_score if `ground_truth_vcf` is not set, and the Pearson correlation otherwise.
Which performance metric should be used can also be specified manually with the `filter_var` parameter.

An example of `snp_filtering_criteria` file is the following:

```
0.8
0.8
0.9
0.9
```

The number of filtering iterations performed is equal to the number of rows in the `snp_filtering_criteria` file.
For each iteration, the filter value of the corresponding line is applied.

## Birneylab-specific information

For ease of use, the ideal settings for stitch for medaka samples have been specified in a profile called `medaka`.
This can be activated with the flag `-profile medaka`.
Always use this profile when working with medaka samples.
For some caveats I cannot set a reference genome directly (medaka in Ensembl is compressed with gzip instead than bgzip) so you still need to provided a path to a valid reference using the `fasta` parameter.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run birneylab/stitchimpute --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> ⚠️ Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run birneylab/stitchimpute -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

<!--
You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).
-->

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull birneylab/stitchimpute
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [birneylab/stitchimpute releases page](https://github.com/birneylab/stitchimpute/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> 💡 If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
