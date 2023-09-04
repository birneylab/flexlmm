include { PLOT_INFO_SCORE } from '../../modules/local/makeplots'
include { PLOT_R2_SITES   } from '../../modules/local/makeplots'
include { PLOT_R2_SAMPLES } from '../../modules/local/makeplots'
include { PLOT_R2_MAF     } from '../../modules/local/makeplots'

workflow PLOTTING {
    take:
    info_score // channel: [mandatory] [ meta, info_score ]
    rsquare    // channel: [optional ] [ meta, r2_per_site, r2_samples, r2_groups ]

    main:
    versions = Channel.empty()
    plots    = Channel.empty()

    rsquare.map { meta, r2_sites, r2_samples, r2_groups -> [ meta, r2_sites   ] }.set { r2_sites   }
    rsquare.map { meta, r2_sites, r2_samples, r2_groups -> [ meta, r2_samples ] }.set { r2_samples }
    rsquare.map { meta, r2_sites, r2_samples, r2_groups -> [ meta, r2_groups  ] }.set { r2_groups  }

    PLOT_INFO_SCORE ( info_score )
    PLOT_R2_SITES   ( r2_sites   )
    PLOT_R2_SAMPLES ( r2_samples )
    PLOT_R2_MAF     ( r2_groups  )

    plots.mix ( PLOT_INFO_SCORE.out.plots ).set { plots }
    plots.mix ( PLOT_R2_SITES.out  .plots ).set { plots }
    plots.mix ( PLOT_R2_SAMPLES.out.plots ).set { plots }
    plots.mix ( PLOT_R2_MAF.out    .plots ).set { plots }

    versions.mix ( PLOT_INFO_SCORE.out.versions ).set { versions }
    versions.mix ( PLOT_R2_SITES  .out.versions ).set { versions }
    versions.mix ( PLOT_R2_SAMPLES.out.versions ).set { versions }
    versions.mix ( PLOT_R2_MAF    .out.versions ).set { versions }

    emit:
    plots    // channel: [ meta, plot ]

    versions // channel: [ versions.yml ]
}
