# Create common plot theme
plot_theme <- ggplot2::theme_bw() +
   ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 20),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.text.y  = ggplot2::element_text(size = 15),
      axis.text.x  = ggplot2::element_text(
         size  = 12,
         angle = 75,
         hjust = 1
      ),
      axis.title.y = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_text(size = 15),
      legend.text  = ggplot2::element_text(size = 15)
   )


#' Make Customer Plots
#'
#' @description Exports a PDF containing a variety of plots meant primarily for
#' distribution to customers. All plots are adapted from their implementation in
#' the Protti R Package (https://cran.r-project.org/web/packages/protti/index.html).
#' By default, coloring for all plots is by condition groups and use the 'Acadia'
#' scheme in the NatParksPalettes R package
#' (https://github.com/kevinsblake/NatParksPalettes).
#'
#' @param sample_correlation A heatmap displaying the spearman rank correlation
#' between each sample. Additionally perform hierarchical clustering between samples.
#' @inheritParams make_qc_plots
#' @param t_test_data The data frame generated from the 'perform_t_test' function.
#' @param volcano A series of volcano plots For each group comparison. log2 fold
#' change are on the x-axis and log10 transformed p-values from the moderated t-
#' test are on the y-axis. By default, fold-change cut-off is set to 1 and p-value
#' significance level is set to 0.05. Any protein groups with values greater than
#' BOTH cutoffs are highlioghted. Plots are faceted by condition.
#' @param fold_change_cutoff The log2 fold change you want to have the volcano plot
#' cut-off at. Default is 1 (Corresponding to a 2x expression change between conditions).
#' @param p_value_cutoff The p-value you want to have the volcano plot cut-off at.
#' Default is 0.05.
#'
#' @return A PDF file written in your working directory
#' @export
#'
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
make_customer_plots <- function(
      pep_data = pep_data,
      pg_data = pg_data,
      t_test_data = t_test_data,
      filename = 'Customer_report.pdf',
      sample_correlation   = T,
      protein_coverage     = T,
      protein_completeness = F,
      protein_intensity    = F,
      volcano              = T,
      fold_change_cutoff   = 1,
      p_value_cutoff       = 0.05
) {

   if(sample_correlation   == T) {

      sample_colors <- 'placeholder'

      ## Calculate the spearman rank correlation between samples --------------
      correlation <- pg_data %>%
         dplyr::distinct(
            sample,
            protein_group,
            log2_pg_maxlfq
         ) %>%
         tidyr::pivot_wider(
            names_from = sample,
            values_from = log2_pg_maxlfq
         ) %>%
         tibble::column_to_rownames(
            var = 'protein_group'
         ) %>%
         stats::cor(
            method = 'spearman',
            use = "complete.obs"
         )

      annotation <- pg_data %>%
         dplyr::mutate(
            condition := as.character(condition)
         ) %>%
         dplyr::distinct(
            sample,
            condition
         ) %>%
         tibble::column_to_rownames(
            var = 'sample'
         )

      ## Create list for coloring of annotations ----------------------------
      n_conditions <- 0
      conditions_colours <- c()

      conditions <- unique(
         dplyr::pull(
            annotation,
            condition
         )
      )

      n_conditions <- length(conditions)

      conditions_colours <- NatParksPalettes::natparks.pals(
         'Denali'
      )[1:n_conditions]

      names(conditions_colours) <- conditions

      annotation_colours <- list(conditions_colours)

      names(annotation_colours) <- 'condition'

      ## Create heatmaply dendrogram for pheatmap -------------------------
      distance <- stats::dist(correlation)
      hierachical_clustering <- stats::hclust(distance)
      dendrogram <- stats::as.dendrogram(hierachical_clustering)
      dendrogram_row <- dendextend::seriate_dendrogram(
         dendrogram,
         distance,
         method = "OLO")
      dendrogram_column <- dendextend::rotate(
         dendrogram_row,
         order = rev(
            labels(distance)[seriation::get_order(
               stats::as.hclust(dendrogram_row)
            )]
         )
      )

      ## Create pheatmap ------------------------------------------------
      sample_heatmap <-
         pheatmap::pheatmap(
            correlation,
            cluster_rows = stats::as.hclust(dendrogram_row),
            cluster_cols = stats::as.hclust(dendrogram_column),
            annotation = annotation,
            annotation_colors = annotation_colours,
            main = "Correlation based hierachical clustering of samples",
            color = NatParksPalettes::natparks.pals('Volcanoes', 100, 'continuous'),
            treeheight_row = 0
         )

   }
   if(protein_coverage     == T) {

      pg_count_data <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            protein_group
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::mutate(
            count = dplyr::n()
         ) %>%
         dplyr::select(-protein_group) %>%
         tidyr::drop_na() %>%
         dplyr::distinct(
            sample,
            condition,
            count
         ) %>%
         dplyr::ungroup()

      pg_count_plot <- pg_count_data %>%
         dplyr::mutate(
            sample := factor(sample,
                             levels = unique(
                                stringr::str_sort(sample,
                                                  numeric = T)
                             )
            )
         ) %>%
         ggplot2::ggplot(
            ggplot2::aes(
               x    = condition,
               y    = count,
               fill = condition,
               group = sample
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1,
            position = 'dodge'
         ) +
         ggplot2::labs(
            title = 'Number of Protein Groups Quantified',
            x     = '',
            y     = 'Count',
            fill  = 'Condition'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               pg_count_data$condition %>%
                  dplyr::n_distinct()
            )
         )
   }
   if(protein_completeness == T) {

      pg_completeness_data <- pep_data %>%
         dplyr::distinct(
            sample,
            protein_group,
            condition,
            log2_pg_maxlfq
         ) %>%
         tidyr::complete(
            sample,
            protein_group
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::summarise(
            completeness =
               sum(!is.na(log2_pg_maxlfq))/dplyr::n() * 100,
            .groups      = 'drop'
         ) %>%
         dplyr::left_join(
            y = pep_data %>%
               dplyr::distinct(
                  sample,
                  condition
               ),
            by = 'sample'
         ) %>%
         dplyr::relocate(condition, .after = sample)

      pg_completeness_plot <- pg_completeness_data %>%
         dplyr::mutate(
            sample := factor(sample,
                             levels = unique(
                                stringr::str_sort(sample,
                                                  numeric = T
                                )
                             )
            )
         ) %>%
         ggplot2::ggplot(
            ggplot2::aes(
               x    = sample,
               y    = completeness,
               fill = condition
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1
         ) +
         ggplot2::geom_text(
            data = pg_completeness_data,
            ggplot2::aes(
               label = round(
                  completeness,
                  digits = 1
               )
            ),
            position = ggplot2::position_stack(vjust = 0.5)
         ) +
         ggplot2::scale_y_continuous(
            limits = c(0, 100)
         ) +
         ggplot2::labs(
            title = 'Protein Group-level Data Completeness',
            x     = '',
            y     = 'Data Completeness [%]'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               pg_completeness_data$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(protein_intensity    == T) {

      pg_intensity_data <- pg_data %>%
         dplyr::distinct(
            sample,
            condition,
            protein_group,
            log2_pg_maxlfq
         ) %>%
         tidyr::drop_na(log2_pg_maxlfq)

      pg_intensity_plot <- pg_intensity_data %>%
         ggplot2::ggplot(
            ggplot2::aes(
               x    = sample,
               y    = log2_pg_maxlfq,
               fill = condition
            )
         ) +
         ggplot2::geom_violin(na.rm = T) +
         ggplot2::geom_boxplot(
            width = 0.15,
            fill  = 'white',
            na.rm = T,
            alpha = 0.6
         ) +
         ggplot2::labs(
            title = 'Protein Group-level Run MaxLFQ Abundances',
            x     = '',
            y     = 'Log2 MaxLFQ'
         ) +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               pg_intensity_data$condition %>%
                  dplyr::n_distinct()
            )
         ) +
         plot_theme

   }
   if(volcano              == T) {

      volcano_data <- t_test_data %>%
         dplyr::mutate(mean_cutoff = p_value_cutoff)

      cutoff_line <- volcano_data %>%
         dplyr::mutate(
            mean_cutoff = -log10(.data$mean_cutoff)
         ) %>%
         dplyr::distinct(
            comparison,
            .data$mean_cutoff
         ) %>%
         tidyr::drop_na(.data$mean_cutoff)

      vol_plot <- volcano_data %>%
         dplyr::filter(
            !(
               (abs(diff) > 1) &
                  ifelse(
                     is.na(pval < .data$mean_cutoff),
                     FALSE,
                     pval < .data$mean_cutoff
                  )
            )
         ) %>%
         ggplot2::ggplot(
            ggplot2::aes(
               label1 = NULL,
               label2 = protein_group
            )
         ) +
         ggplot2::geom_point(
            ggplot2::aes(
               x = diff,
               y = -log10(pval)
            ),
            colour = '#9da7bf'
         ) +
         ggplot2::geom_point(
            data = dplyr::filter(
               volcano_data,
               (abs(diff) > fold_change_cutoff) &
                  (pval       < .data$mean_cutoff)
            ),
            ggplot2::aes(
               x = diff,
               y = -log10(pval)
            ),
            size   = 3,
            colour = '#006375'
         ) +
         ggplot2::labs(
            title        = 'Volcano Plot',
            x_axis_label = 'log2(fold change)',
            y_axis_label = '-log10(p-value)'
         ) +
         ggplot2::geom_hline(
            data = cutoff_line,
            ggplot2::aes(
               yintercept = .data$mean_cutoff
            ),
            linetype = 'dashed'
         ) +
         ggplot2::geom_vline(
            xintercept = fold_change_cutoff,
            linetype   = 'dashed'
         ) +
         ggplot2::geom_vline(
            xintercept = -fold_change_cutoff,
            linetype   = 'dashed'
         ) +
         ggplot2::theme_bw() +
         ggplot2::theme(
            plot.title       = ggplot2::element_text(size = 20),
            axis.title.x     = ggplot2::element_text(size = 15),
            axis.text.y      = ggplot2::element_text(size = 15),
            axis.text.x      = ggplot2::element_text(size = 12),
            axis.title.y     = ggplot2::element_text(size = 15),
            strip.text       = ggplot2::element_text(size = 20),
            strip.background = ggplot2::element_blank(),
            legend.position  = 'none'
         ) +
         ggplot2::scale_x_continuous(
            breaks = seq(
               round(
                  -1 * max(abs(volcano_data %>%
                                  dplyr::select(diff)), na.rm = T) - 0.5, 0
               ),
               round(
                  max(abs(volcano_data %>%
                             dplyr::select(diff)), na.rm = T) + 0.5, 0
               ),
               1
            )
         ) +
         ggplot2::coord_cartesian(
            xlim = c(
               round(
                  -1 * max(abs(volcano_data %>%
                                  dplyr::select(diff)), na.rm = T) - 0.5, 0
               ),
               round(
                  max(abs(volcano_data %>%
                             dplyr::select(diff)), na.rm = T) + 0.5, 0
               )
            )
         ) +
         ggforce::facet_wrap_paginate(~ comparison, ncol = 1, nrow = 1)

   }

   # Print the created plots
   pdf(file = filename)
   if(sample_correlation     == T) {
      print(sample_heatmap)
   }
   if(protein_coverage       == T) {
      print(pg_count_plot)
   }
   if(protein_completeness   == T) {
      print(pg_completeness_plot)
   }
   if(protein_intensity      == T) {
      print(pg_intensity_plot)
   }
   if(volcano                == T) {
      for(i in 1:n_pages(vol_plot)) {
         vol_save <- vol_plot +
            facet_wrap_paginate(~ comparison, ncol = 1, nrow = 1, page = i)
         print(vol_save)
      }
   }

}
