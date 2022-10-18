pep_count_fun <- function(df) {
   as.data.frame(
      table(
         unique(
            df[,
               c(
                  'protein_group',
                  'modified_peptide'
               )
            ]
         )$'protein_group'
      )
   )
}

#' Combine Metadata with Peptide-level Data
#'
#' @description Annotates peptide-level data with metadata and counts the number of times a peptide was quantified
#'
#' @param data The data frame output from the prepare_data function
#' @param sample_annotation The matrix output from the create_metadata function
#'
#' @return A data frame containing peptide-level MaxLFQ intensities and peptide count for each protein group
#' @export
#'
annotate_peptide <- function(
      data,
      sample_annotation) {

   pep_data_filtered <- data %>%
      dplyr::left_join(
         y = sample_annotation,
         by = 'sample'
      ) %>%
      dplyr::relocate(
         condition,
         .after = sample
      )

   #Count the number of peptides quantified (per protein group) in each samples
   pep_count <- data %>%
      dplyr::select(
         c(
            'sample',
            'protein_group',
            'modified_peptide'
         )
      ) %>%
      dplyr::group_by(sample) %>%
      tidyr::nest() %>%
      dplyr::mutate(
         n_peptides = purrr::map(data, pep_count_fun)
      ) %>%
      dplyr::select(
         -c('data')
      ) %>%
      tidyr::unnest(n_peptides)  %>%
      purrr::set_names(
         'sample',
         'protein_group',
         'number_precursors'
      )

   pep_count$sample <- sub("^[^_]*_([^_]*).*", "\\1", pep_count$sample)

   peptide_data <- pep_data_filtered %>%
      dplyr::left_join(
         y = pep_count,
         by = c('sample', 'protein_group')
      ) %>%
      dplyr::relocate(
         number_precursors,
         .after = peak_width
      )

   return(peptide_data)
}
