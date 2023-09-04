// rename a file according to prefix, maintaining the same extension
process RENAME_VCF_TBI {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path(new_vcf), path(new_tbi) , emit: renamed
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args   ?: ""
    def vcf_ext    = ( vcf.name - vcf.simpleName )
    def tbi_ext    = ( tbi.name - tbi.simpleName )
    new_vcf        = prefix + vcf_ext
    new_tbi        = prefix + tbi_ext
    """
    ln -sf \$PWD/${vcf} ${new_vcf}
    ln -sf \$PWD/${tbi} ${new_tbi}

    touch versions.yml
    """
}
