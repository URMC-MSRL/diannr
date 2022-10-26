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

impute_normal <- function(
      df,
      width     = 0.3,
      downshift = 1.8,
      seed      = 100) {

   if (!is.matrix(df))
      df <- as.matrix(df)
   mx <- max(df, na.rm = T)
   mn <- min(df, na.rm = T)
   if (mx - mn > 20)
      warning("Please make sure the values are log-transformed")

   set.seed(seed)
   df <- apply(
      df,
      2,
      function(temp) {
         temp[!is.finite(temp)] <- NA
         temp_sd <- stats::sd(temp, na.rm = T)
         temp_mean <- mean(temp, na.rm    = T)
         shrinked_sd <- width * temp_sd
         downshifted_mean <- temp_mean - downshift * temp_sd
         n_missing <- sum(is.na(temp))
         temp[is.na(temp)] <- stats::rnorm(
            n_missing,
            mean = downshifted_mean,
            sd   = shrinked_sd)
         temp
      })
   return(df)
}

#' Combine Metadata with Protein-level Data
#'
#' @description Annotates protein group-level data with metadata, imputes
#' missing values using a standard distribution as implemented in Perseus, and
#' counts the number of times a peptide was quantified
#'
#' @inheritParams annotate_peptide
#'
#' @return A data frame containing protein-level MaxLFQ intensities and peptide count for each protein group
#' @export
#'
annotate_protein <- function(
      data = data,
      sample_annotation = sample_annotation
) {
   ### Count the number of peptides quantified ---------------------------------
   pep_count <- data %>%
      dplyr::select(
         c('sample',
           'protein_group',
           'modified_peptide')
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

   ### Impute missing values ----------------------------------------------------
   imputed_values <- data %>%
      dplyr::distinct(
         sample,
         protein_group,
         log2_pg_maxlfq
      ) %>%
      tidyr::pivot_wider(
         names_from = sample,
         values_from = log2_pg_maxlfq
      ) %>%
      tibble::column_to_rownames(var = 'protein_group') %>%
      impute_normal(.)

   imputed_with_names <- dplyr::bind_cols(
      imputed_values %>%
         as.data.frame() %>%
         tibble::rownames_to_column(var = 'protein_group'),
      data %>%
         dplyr::distinct(
            sample,
            protein_group,
            gene_name,
            protein_name,
            log2_pg_maxlfq
         ) %>%
         tidyr::pivot_wider(
            names_from = sample,
            values_from = log2_pg_maxlfq
         ) %>%
         dplyr::select(
            c(gene_name,
              protein_name)
         )
   ) %>%
      dplyr::relocate(
         c(gene_name:protein_name),
         .after = protein_group
      ) %>%
      tidyr::pivot_longer(
         cols = 4:(3 + dplyr::n_distinct(sample_annotation)),
         names_to = 'sample',
         values_to = 'log2_pg_maxlfq'
      ) %>%
      dplyr::relocate(
         sample,
         .before = protein_group
      ) %>%
      dplyr::left_join(
         y = sample_annotation,
         by = 'sample'
      ) %>%
      dplyr::relocate(
         condition,
         .after = sample
      )

   ### Combine imputed quantities with peptide count ----------------------------
   protein_data <- imputed_with_names %>%
      dplyr::left_join(
         y = pep_count,
         by = c('sample', 'protein_group'),
         na_matched = 'na'
      )

   return(protein_data)

}
