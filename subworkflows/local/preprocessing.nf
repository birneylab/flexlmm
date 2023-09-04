//
// Create cramlist for STITCH, index reference fasta, and get chromosome names from fasta
// index
//

include { DOWNSAMPLE                                    } from '../../subworkflows/local/downsample'
include { SAMTOOLS_FAIDX                                } from '../../modules/nf-core/samtools/faidx'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_GROUND_TRUTH } from '../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FREQ_VCF     } from '../../modules/nf-core/bcftools/index'

workflow PREPROCESSING {
    take:
    reads               // channel: [mandatory] [ meta, cram, crai ]
    fasta               // channel: [mandatory] [ fasta ]
    skip_chr            // channel: [mandatory] list of chromosomes to skip
    ground_truth        // channel: [optional]  [ vcf_file ]
    freq_vcf            // channel: [optional]  [ vcf_file ]
    downsample_coverage // channel: [optional]  coverage to downsample ground truth to

    main:
    versions = Channel.empty()

    fasta.map { fasta -> [ [ id:fasta.baseName ], fasta ] }.set { fasta }
    SAMTOOLS_FAIDX ( fasta, [['id':null], []] )

    fasta
    .join ( SAMTOOLS_FAIDX.out.fai, failOnMismatch: true, failOnDuplicate: true )
    .first ()
    .set { reference }

    reference
    .map{ meta, fasta, fasta_fai -> fasta_fai }
    .splitCsv ( sep:"\t" )
    .map { name, length, offset, linebases, linewidth -> name }
    .filter { name -> ! skip_chr.contains ( name ) }
    .collect ()
    .ifEmpty ( ["stub_chr"] )
    .set { chr_list }

    if ( params.downsample_coverage ) {
        DOWNSAMPLE( reads, skip_chr, downsample_coverage )
        DOWNSAMPLE.out.reads.set { reads }

        versions.mix ( DOWNSAMPLE.out.versions ) .set { versions }
    }

    reads
    .map { meta, cram, crai -> cram[-1] as String } // cram name without path
    .collectFile ( name: "stitch_cramlist.txt", newLine:true )
    .set { stitch_cramlist }

    reads
    .map { meta, cram, crai -> [["id": "collected_samples"], cram, crai] }
    .groupTuple ()
    .combine ( stitch_cramlist )
    .collect () // needed to make it broadcastable
    .set { collected_samples }

    ground_truth.map { [["id": "ground_truth_vcf"], it] }.set { ground_truth }
    BCFTOOLS_INDEX_GROUND_TRUTH ( ground_truth )

    ground_truth
    .join( BCFTOOLS_INDEX_GROUND_TRUTH.out.csi, failOnMismatch: true, failOnDuplicate: true )
    .set { ground_truth }

    freq_vcf.map { [["id": "freq_vcf"], it] }.set { freq_vcf }
    BCFTOOLS_INDEX_FREQ_VCF ( freq_vcf )

    freq_vcf
    .join( BCFTOOLS_INDEX_FREQ_VCF.out.csi, failOnMismatch: true, failOnDuplicate: true )
    .set { freq_vcf }

    // Gather versions of all tools used
    versions.mix ( SAMTOOLS_FAIDX.out.versions              ) .set { versions }
    versions.mix ( BCFTOOLS_INDEX_GROUND_TRUTH.out.versions ) .set { versions }

    emit:
    collected_samples // channel: [ meta, collected_crams, collected_crais, stitch_cramlist ]
    reference         // channel: [ meta, fasta, fasta_fai ]
    chr_list          // channel: list of chromosomes to consider
    ground_truth      // channel: [ meta, vcf, vcf_index ]
    freq_vcf          // channel: [ meta, vcf, vcf_index ]

    versions          // channel: [ versions.yml ]
}
