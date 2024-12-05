process BGEN_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(bgen), path(gen), path(sample)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    def use_dosage_flags = use_dosage ? "" : "fill-missing-from-dosage erase-dosage"
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --bgen ${bgen} --gen ${gen} --sample ${sample} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

process BCF_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(bcf)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --bcf ${bcf} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

process VCF_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(vcf)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    def use_dosage_flags = use_dosage ? "" : "fill-missing-from-dosage erase-dosage"
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --vcf ${vcf} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

process PGEN_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    def use_dosage_flags = use_dosage ? "" : "fill-missing-from-dosage erase-dosage"
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --pgen ${pgen} --psam ${psam} --pvar ${pvar} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

process BED_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    def use_dosage_flags = use_dosage ? "" : "fill-missing-from-dosage erase-dosage"
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --bed ${bed} --bim ${bim} --fam ${fam} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}

process PED_TO_PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a3.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a3.7--h4ac6f70_4' :
        'biocontainers/plink2:2.00a3.7--h4ac6f70_4' }"

    input:
    tuple val(meta), path(ped), path(map)
    val maf_min
    val use_dosage

    output:
    tuple val(meta), path("*.pgen"), path("*.pvar.zst"), path("*.psam"), emit: pgen_psam_pvar
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def args3      = task.ext.args3 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def mem_mb     = task.memory.toMega()
    def maf_filter = maf_min ? "--maf ${maf_min}" : ""
    def use_dosage_flags = use_dosage ? "" : "fill-missing-from-dosage erase-dosage"
    """
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --set-all-var-ids @_#_\\\$r_\\\$a \\
        --min-alleles 2 \\
        --max-alleles 2 \\
        $maf_filter \\
        --out $prefix \\
        $args \\
        --ped ${ped} --map ${map} ${args2} \\
        --make-pgen vzs $args3 $use_dosage_flags

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_mb = task.memory.toMega()
    """
    touch ${prefix}.pgen
    touch ${prefix}.psam
    touch ${prefix}.pvar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}




