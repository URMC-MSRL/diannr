# Calculate group mean --------------------------------------------------------
mean_function <- function(df) {
   df %>%
      dplyr::group_by(protein_group) %>%
      dplyr::mutate(mean_log2_pg_maxlfq = mean(log2_pg_maxlfq)) %>%
      dplyr::distinct(protein_group,
                      mean_log2_pg_maxlfq)
}

calculate_group_means <- function(df) {
   mean_log2_data <- df %>%
      dplyr::select(c(condition,
                      protein_group,
                      log2_pg_maxlfq)) %>%
      dplyr::group_by(condition) %>%
      tidyr::nest() %>%
      dplyr::mutate(mean_data = purrr::map(data,
                                           mean_function),
                    .keep = 'unused')  %>%
      tidyr::unnest(col = mean_data)

   return(mean_log2_data)
}

#' Generate Customer Report
#'
#' @description Exports a .tsv into your working directory that can be further
#' edited and sent as a report to customers. Each row corresponds to a different
#' protein group and has rows for:
#'     MaxLFQ protein abundance of each sample (in log2)
#'     Number of peptides quantified for each sample
#'     Average MaxLFQ for each group
#'     Log2 fold-change of each group comparison
#'     p-value of each group comparison
#'
#' @param t_test_data The data frame output from the 'perform_t_test' function.
#' @param pg_data The data frame output from the 'annotate_protein' function.
#' @param filename The desired name of the exported report. Should be in quotes
#' and end in .tsv
#'
#' @return A .tsv file written in your working directory.
#' @export
#'
#' @importFrom readr write_tsv
#'
write_report <- function(t_test_data = t_test_data,
                         pg_data = pg_data,
                          filename = 'output_report.tsv') {

   group_means <- calculate_group_means(pg_data) %>%
      tidyr::pivot_wider(names_from = 'condition',
                         values_from = 'mean_log2_pg_maxlfq',
                         names_prefix = 'mean_maxlfq_')

   report_data <- dplyr::full_join(x = dplyr::full_join(t_test_data %>%
         dplyr::select(-c(`CI_2.5`:avg_abundance)) %>%
         dplyr::relocate(comparison,
                         .after = protein_group) %>%
         dplyr::rename('log2_FC' = diff,
                       'p_value' = pval) %>%
         tidyr::pivot_wider(names_from = comparison,
                            values_from = c(log2_FC:p_value)),
                                                        pg_data %>%
         dplyr::distinct(sample,
                         condition,
                         protein_group,
                         log2_pg_maxlfq,
                         number_precursors) %>%
         dplyr::select(-condition) %>%
         dplyr::rename(number_peptides = number_precursors) %>%
         tidyr::pivot_wider(names_from = sample,
                            values_from = c(number_peptides:log2_pg_maxlfq)),
                                                        by = 'protein_group'),
                                    y = group_means %>%
         dplyr::left_join(x = .,
                          y = pg_data %>%
                             dplyr::distinct(protein_group,
                                             protein_name,
                                             gene_name),
                          by = 'protein_group')) %>%
      dplyr::relocate(c(protein_name,
                        gene_name),
                      .after = protein_group)
   readr::write_tsv(report_data, 
                    filename)
}
