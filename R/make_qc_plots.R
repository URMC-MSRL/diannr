# Create common plote theme
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

#' Make QC Plots
#'
#' @description Exports a PDF containing a variety of plots meant primarily for
#' internal QC purposes to quickly assess quality of the data. All plots are adapted
#' from their implementation in the Protti R Package
#' (https://cran.r-project.org/web/packages/protti/index.html). By default,
#' coloring for all plots is by condition groups and use the 'Acadia' scheme in the
#' NatParksPalettes R package (https://github.com/kevinsblake/NatParksPalettes).
#'
#' @param pep_data The data frame generated from the 'annotate_peptide' function.
#' @param pg_data The data frame generated from the 'annotate_protein' function
#' @param filename The desired name of the exported PDF. Should be in quotes and
#' end in .pdf
#' @param peak_shape A line plot depicting the mean peak width for each sample as
#' a function of retention time.
#' @param charge_state A bar graph depicting the percentage of precursor charge
#' states for every sample.
#' @param cv A box and whisker plot over a violin plot depicting the distribution
#' of protein-group MaxLFQ abundances among:
#'     all samples together
#'     all samples in each conditions
#' @param precursor_coverage A bar graph depicting the number of precursors
#' quantified in a given sample.
#' @param peptide_coverage A bar graph depicting the number of peptides (can
#' include multiple charge states) quantified in a given sample.
#' @param protein_coverage A bar graph depicting the number of proteins
#' quantified in a given sample.
#' @param precursor_completeness A bar graph depicting the percentage of precursors
#' quantified in a given sample as compared to all precursors quantified.
#' @param peptide_completeness  A bar graph depicting the percentage of peptides
#' (can include multiple charge states) quantified in a given sample as compared
#' to all peptides quantified.
#' @param protein_completeness A bar graph depicting the percentage of proteins
#' quantified in a given sample as compared to all proteins quantified.
#' @param precursor_intensity A box and whisker plot over a violin plot depicting
#' the distribution of normalized precursor intensities of each sample.
#' @param protein_intensity A box and whisker plot over a violin plot depicting
#' the distribution of MaxLFQ protein abundances of each sample.
#'
#' @return A PDF file written in your working directory.
#' @export
#'
make_qc_plots <- function(
      pep_data = pep_data,
      pg_data = pg_data,
      filename = 'QC_report.pdf',
      peak_shape             = T,
      charge_state           = T,
      cv                     = T,
      precursor_coverage     = T,
      peptide_coverage       = T,
      protein_coverage       = T,
      precursor_completeness = T,
      peptide_completeness   = T,
      protein_completeness   = T,
      precursor_intensity    = T,
      protein_intensity      = T
) {

   if(peak_shape             == T) {

      peak_width_data <- pep_data %>%
         dplyr::distinct(
            sample,
            rt,
            peak_width
         )

      peak_width_plot <- peak_width_data %>%
         dplyr::mutate(
            sample := factor(sample,
                             levels = unique(
                                stringr::str_sort(
                                   sample,
                                   numeric = T
                                )
                             )
            )
         ) %>%
         ggplot2::ggplot(
            ggplot2::aes(
               x = rt,
               y = peak_width
            )
         ) +
         ggplot2::stat_summary_bin(
            ggplot2::aes(
               col = sample
            ),
            size     = 1,
            geom     = 'line',
            binwidth = 1,
            fun      = median
         ) +
         ggplot2::labs(
            title = 'Median Peak Width Over Retention Time',
            x     = 'Retention Time [min]',
            y     = 'Median Peak Width [min]',
            color = 'Sample'
         ) +
         ggplot2::scale_color_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               peak_width_data$sample %>%
                  dplyr::n_distinct()
            )
         ) +
         plot_theme


   }
   if(charge_state           == T) {

      charge_state_data <- pep_data %>%
         tidyr::drop_na(precursor_normalized) %>%
         dplyr::distinct(
            sample,
            precursor,
            charge
         ) %>%
         dplyr::count(
            sample,
            charge
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::mutate(
            total_peptides = sum(n)
         ) %>%
         dplyr::group_by(
            sample,
            charge
         ) %>%
         dplyr::summarise(
            charge_per = n/total_peptides * 100
         ) %>%
         dplyr::ungroup() %>%
         dplyr::mutate(
            charge := forcats::fct_inorder(
               factor(charge)
            )
         ) %>%
         dplyr::mutate(
            sample := factor(
               sample,
               levels = unique(
                  stringr::str_sort(
                     sample,
                     numeric = T
                  )
               )
            )
         )

      charge_state_plot <- charge_state_data %>%
         ggplot2::ggplot(
            aes(
               x    = sample,
               y    = charge_per,
               fill = charge,
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1
         ) +
         ggplot2::geom_text(
            data = charge_state_data %>%
               dplyr::filter(
                  charge_per < 5
               ),
            ggplot2::aes(
               label = round(charge_per,
                             digits = 1
               )
            ),
            position = ggplot2::position_stack(vjust = 0.9)
         ) +
         ggplot2::labs(
            title    = 'Charge Distribution',
            subtitle = 'By Percent of Total Peptide Count',
            x        = '',
            y        = '% of Total Peptide Count',
            fill     = 'Charge'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               charge_state_data$charge %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(cv                     == T) {

      cv_data <- pg_data %>%
         dplyr::distinct(
            protein_group,
            condition,
            log2_pg_maxlfq
         ) %>%
         dplyr::rename(
            'intensity' = log2_pg_maxlfq
         ) %>%
         dplyr::mutate(
            'intensity' = 2^intensity,
            .keep       = 'unused'
         ) %>%
         tidyr::drop_na(intensity) %>%
         dplyr::group_by(protein_group) %>%
         dplyr::mutate(
            cv_combined = (stats::sd(intensity)/mean(intensity)) * 100
         ) %>%
         dplyr::group_by(
            condition, protein_group
         ) %>%
         dplyr::mutate(
            cv = (stats::sd(intensity)/mean(intensity)) * 100
         ) %>%
         dplyr::ungroup() %>%
         dplyr::distinct(
            condition,
            protein_group,
            cv_combined,
            cv
         ) %>%
         drop_na() %>%
         tidyr::pivot_longer(
            cols = starts_with('cv'),
            names_to  = 'type',
            values_to = 'values'
         ) %>%
         dplyr::mutate(
            type = ifelse(
               type == 'cv', condition, 'combined'
            )
         ) %>%
         dplyr::mutate(
            type = forcats::fct_relevel(
               as.factor(type),
               'combined')
         ) %>%
         dplyr::select(-condition) %>%
         dplyr::group_by(type) %>%
         dplyr::mutate(
            median = stats::median(values)
         ) %>%
         dplyr::distinct() %>%
         dplyr::filter(values < 200)

      cv_plot <- ggplot2::ggplot(
         cv_data,
         ggplot2::aes(
            x    = type,
            y    = values,
            fill = type
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
            title = 'Protein Group Coefficients of Variation',
            x     = '',
            y     = 'Coefficient of variation [%]',
            fill  = 'Condition'
         ) +
         ggplot2::scale_y_continuous(
            limits = c(0, 200)
         ) +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               cv_data$type %>%
                  dplyr::n_distinct()
            )
         ) +
         plot_theme


   }
   if(precursor_coverage     == T) {

      precursor_count <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            precursor
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::mutate(
            count = dplyr::n()
         ) %>%
         dplyr::select(-c('precursor')) %>%
         dplyr::distinct() %>%
         dplyr::ungroup()

      precursor_count_plot <- precursor_count %>%
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
               x    = sample,
               y    = count,
               fill = condition
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1
         ) +
         ggplot2::labs(
            title = 'Number of Precursors Quantified',
            x     = '',
            y     = 'Count',
            fill  = 'Condition'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               precursor_count$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(peptide_coverage       == T) {

      peptide_count <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            peptide
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::mutate(
            count = dplyr::n()
         ) %>%
         dplyr::select(-peptide) %>%
         dplyr::distinct() %>%
         dplyr::ungroup()

      peptide_count_plot <- peptide_count%>%
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
               x    = sample,
               y    = count,
               fill = condition
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1
         ) +
         ggplot2::labs(
            title = 'Number of Peptides Quantified',
            x     = '',
            y     = 'Count',
            fill  = 'Condition'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               peptide_count$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(protein_coverage       == T) {

      pg_count <- pep_data %>%
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

      pg_count_plot <- pg_count %>%
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
               x    = sample,
               y    = count,
               fill = condition
            )
         ) +
         ggplot2::geom_col(
            col  = 'black',
            size = 1
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
               pg_count$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(precursor_completeness == T) {

      precursor_completeness_data <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            precursor,
            log2_precursor_normalized
         ) %>%
         tidyr::complete(
            sample,
            precursor
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::summarise(
            completeness =
               sum(!is.na(log2_precursor_normalized))/dplyr::n() * 100,
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

      precursor_completeness_plot <- precursor_completeness_data %>%
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
            fill = NatParksPalettes::natparks.pals('Denali', 1),
            col  = 'black',
            size = 1
         ) +
         ggplot2::geom_text(
            data = precursor_completeness_data,
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
            title = 'Precursor-level Data Completeness',
            x     = '',
            y     = 'Data Completeness [%]'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               precursor_completeness_data$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(peptide_completeness   == T) {

      peptide_completeness_data <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            peptide,
            log2_precursor_normalized
         ) %>%
         tidyr::complete(
            sample,
            peptide
         ) %>%
         dplyr::group_by(sample) %>%
         dplyr::summarise(
            completeness =
               sum(!is.na(log2_precursor_normalized))/dplyr::n() * 100,
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

      peptide_completeness_plot <- peptide_completeness_data %>%
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
            data = peptide_completeness_data,
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
            title = 'Peptide-level Data Completeness',
            x     = '',
            y     = 'Data Completeness [%]'
         ) +
         plot_theme +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               peptide_completeness_data$condition %>%
                  dplyr::n_distinct()
            )
         )

   }
   if(protein_completeness   == T) {

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
         )  %>%
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
   if(precursor_intensity    == T) {

      precursor_intensity_data <- pep_data %>%
         dplyr::distinct(
            sample,
            condition,
            precursor,
            log2_precursor_normalized
         ) %>%
         tidyr::drop_na(log2_precursor_normalized)

      precursor_intensity_plot <- precursor_intensity_data %>%
         ggplot2::ggplot(
            ggplot2::aes(
               x    = sample,
               y    = log2_precursor_normalized,
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
            title = 'Precursor-level Run Normalized Intensities',
            x     = '',
            y     = 'Log2 Normalized Intensity'
         ) +
         ggplot2::scale_fill_manual(
            values = NatParksPalettes::natparks.pals(
               'Denali',
               precursor_intensity_data$condition %>%
                  dplyr::n_distinct()
            )
         ) +
         plot_theme

   }
   if(protein_intensity      == T) {

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

   pdf(file = filename)
   if(peak_shape             == T) {
      print(peak_width_plot)
   }
   if(charge_state           == T) {
      print(charge_state_plot)
   }
   if(cv                     == T) {
      print(cv_plot)
   }
   if(precursor_coverage     == T) {
      print(precursor_count_plot)
   }
   if(peptide_coverage       == T) {
      print(peptide_count_plot)
   }
   if(protein_coverage       == T) {
      print(pg_count_plot)
   }
   if(precursor_completeness == T) {
      print(precursor_completeness_plot)
   }
   if(peptide_completeness   == T) {
      print(peptide_completeness_plot)
   }
   if(protein_completeness   == T) {
      print(pg_completeness_plot)
   }
   if(precursor_intensity    == T) {
      print(precursor_intensity_plot)
   }
   if(protein_intensity      == T) {
      print(pg_intensity_plot)
   }
   dev.off()
}
