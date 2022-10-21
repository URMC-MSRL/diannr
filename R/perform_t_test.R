#' Perform Moderated T-Test
#'
#' @description Calculates differential abundance and performs a moderated
#' t-test as implemented in the protti R package (https://cran.r-project.org/web/packages/protti/index.html).
#' Calculations are performed on protein-group abundances and between all possible
#' comparisons.
#' @param pg_data Data frame output from the 'annotate_protein' function
#'
#' @return A data frame containing the log2 fold change between conditions for each protein group and its associated p-value.
#' @export
#'
#' @importFrom rlang :=
#'
perform_t_test <- function(
      pg_data = pg_data
      ) {
   ## Assign Missing-ness ------------------------------------------------------
   message('[1/4] Assigning missing-ness ...', appendLF = F)
   ### Create all pair-wise comparisons -----------------------------------------
   all_conditions <- unique(
      dplyr::pull(
         pg_data,
         condition
      )
   )

   all_combinations <- tibble::as_tibble(
      t(
         utils::combn(
            all_conditions,
            m = 2
         )
      ),
      .name_repair = ~ make.names(.,
                                  unique = T)
   ) %>%
      dplyr::rename(
         V1 = .data$X,
         V2 = .data$X.1
      ) %>%
      dplyr::mutate(
         combinations = paste0(
            .data$V1,
            '_vs_',
            .data$V2
         )
      ) %>% #Create a data frame that contains all combinations to be tested
      tidyr::pivot_longer(
         cols = c(
            .data$V1,
            .data$V2
         ),
         names_to = 'name',
         values_to = 'condition'
      ) %>%
      dplyr::select(-.data$name) %>%
      dplyr::group_by(condition) %>%
      dplyr::mutate(
         comparison = list(.data$combinations)
      ) %>%
      dplyr::distinct(
         .data$comparison,
         condition
      )

   data_prep <- data %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
         sample,
         condition,
         protein_group,
         log2_pg_maxlfq
      ) %>%
      tidyr::complete(
         nesting(
            sample,
            condition
         ),
         protein_group
      ) %>%
      dplyr::group_by(
         protein_group,
         condition
      ) %>%
      dplyr::mutate(
         n_detect = sum(!is.na(log2_pg_maxlfq))
      ) %>%
      dplyr::mutate(
         n_replicates = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
         all_combinations,
         by = 'condition'
      ) %>%
      tidyr::unnest(.data$comparison)

   ### Unequal replicate comparison check --------------------------------------
   unequal_replicates <- data_prep %>%
      dplyr::arrange(condition) %>%
      dplyr::distinct(
         .data$n_replicates,
         .data$comparison
      ) %>%
      dplyr::group_by(.data$comparison) %>%
      dplyr::mutate(
         n = dplyr::n()
      ) %>%
      dplyr::filter(.data$n > 1) %>%
      dplyr::mutate(
         n_replicates = paste0(
            .data$n_replicates,
            collapse = '/'
         )
      ) %>%
      dplyr::distinct(
         .data$n_replicates,
         .data$comparison
      )

   if (
      nrow(unequal_replicates) != 0) {
      message("\n")
      message(
         strwrap(
            "The following comparisons have been detected to have unequal
          replicate numbers. If this is intended please ignore this message.
          This function can appropriately deal with unequal replicate numbers.",
          prefix = "\n",
          initial = ""),
         "\n"
      )
      message(
         paste0(
            utils::capture.output(
               unequal_replicates
            ),
            collapse = "\n"
         )
      )
   }

   ### Assign missing values ---------------------------------------------------
   comparison_data <- data_prep %>%
      dplyr::mutate(
         type = ifelse(
            condition == stringr::str_extract(
               .data$comparison,
               pattern = '(?<=_vs_).+'),
            'control',
            'experimental'
         )
      ) %>%
      split(.$comparison) %>%
      purrr::map_df(
         .f = ~ .x %>%
            tidyr::pivot_wider(
               names_from = .data$type,
               values_from = c(
                  .data$n_detect,
                  .data$n_replicates
               )
            ) %>%
            dplyr::group_by(protein_group) %>%
            tidyr::fill(
               .data$n_detect_experimental,
               .data$n_detect_control,
               .data$n_replicates_experimental,
               .data$n_replicates_control,
               .direction = 'updown'
            ) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(
               missingness = dplyr::case_when(
                  .data$n_detect_control == .data$n_replicates_control &
                     .data$n_detect_experimental == .data$n_replicates_experimental ~ 'complete',
                  .data$n_detect_control <= floor(
                     .data$n_replicates_control * 0.2
                  ) &
                     .data$n_detect_experimental == .data$n_replicates_experimental ~ 'MNAR',
                  .data$n_detect_control == .data$n_replicates_control &
                     .data$n_detect_experimental <= floor(
                        n_replicates_experimental * 0.2
                     ) ~ 'MNAR',
                  .data$n_detect_control >= max(
                     floor(.data$n_replicates_control * 0.7), 1
                  ) &
                     .data$n_detect_experimental >= max(
                        floor(.data$n_replicates_control * 0.7), 1
                     ) ~ 'MAR'
               )
            )
      ) %>%
      dplyr::select(
         -c(
            .data$n_detect_control,
            .data$n_detect_experimental,
            .data$n_replicates_control,
            .data$n_replicates_experimental
         )
      ) %>%
      dplyr::arrange(
         factor(
            protein_group,
            levels = unique(
               stringr::str_sort(protein_group),
               numeric = T
            )
         )
      )
   message('DONE', appendLF = F)

   ## t-test -------------------------------------------------------------------
   ### Create t-test input data ------------------------------------------------
   message('[2/4] Defining the parameters for the moderated t-test ...',
           appendLF = F)
   moderated_t_test_input <- comparison_data %>%
      dplyr::distinct(
         protein_group,
         sample,
         log2_pg_maxlfq
      ) %>%
      tidyr::drop_na(log2_pg_maxlfq) %>%
      tidyr::pivot_wider(
         names_from = sample,
         values_from = log2_pg_maxlfq
      ) %>%
      tibble::column_to_rownames(
         var = rlang::as_name('protein_group')
      ) %>%
      as.matrix()

   ### Defining the t-test design ----------------------------------------------
   moderated_t_test_map <- comparison_data %>%
      dplyr::distinct(
         sample,
         condition
      ) %>%
      dplyr::mutate(
         condition := paste0('x', condition)
      ) %>%
      dplyr::arrange(sample)

   moderated_t_test_design <- stats::model.matrix(
      ~ 0 + factor(
         stringr::str_replace_all(
            dplyr::pull(
               moderated_t_test_map,
               condition
            ),
            pattern = ' ',
            replacement = '_'
         )
      )
   )

   colnames(moderated_t_test_design) <- levels(
      factor(
         stringr::str_replace_all(
            dplyr::pull(
               moderated_t_test_map,
               condition
            ),
            pattern = ' ',
            replacement = '_'
         )
      )
   )
   message('DONE', appendLF = F)

   ### Fitting the regression model --------------------------------------------
   message('[3/4] Performing the moderated t-test', appendLF = F)

   moderated_t_test_fit <- suppressWarnings(
      limma::lmFit(
         moderated_t_test_input,
         moderated_t_test_design)
   )

   ### Create a matrix of differences ------------------------------------------
   names <- paste0(
      'x',
      stringr::str_extract(
         unique(
            dplyr::pull(
               comparison_data,
               comparison)),
         pattern = '.+(?=_vs_)'), '_vs_x',
      stringr::str_extract(
         unique(
            dplyr::pull(
               comparison_data, comparison
            )
         ),
         pattern = '(?<=_vs_).+'
      )
   )

   comparisons <- paste0(
      'x',
      stringr::str_extract(
         stringr::str_replace_all(
            unique(
               dplyr::pull(
                  comparison_data,
                  comparison
               )
            ),
            pattern = ' ',
            replacement = '_'
         ),
         pattern = '.+(?=_vs_)'
      ),
      '-x',
      stringr::str_extract(
         stringr::str_replace_all(
            unique(
               dplyr::pull(
                  comparison_data,
                  comparison
               )
            ),
            pattern = ' ',
            replacement = '_'
         ),
         pattern = '(?<=_vs_).+'
      )
   )

   combinations <- purrr::map2(
      .x = names,
      .y = comparisons,
      .f = function(x, y) {
         rlang::exprs(
            !!rlang::as_name(x) := !!y
         )
      }
   )

   contrast_matrix <- eval(
      rlang::expr(
         limma::makeContrasts(
            !!!unlist(combinations),
            levels = moderated_t_test_design
         )
      )
   )

   ### Compute Contrasts From Regression Model ---------------------------------
   moderated_t_test_fit2 <- limma::contrasts.fit(
      moderated_t_test_fit,
      contrast_matrix)

   ### Compute Empirical Bayes Statistics --------------------------------------
   moderated_t_test_fit3 <- limma::eBayes(
      moderated_t_test_fit2
   )
   message('DONE', appendLF = F)

   ### Create Result Table -----------------------------------------------------
   message('[4/4] Creating a data frame from the t-test results ...',
           appendLF = F)

   moderated_test_missingness <- comparison_data %>%
      tidyr::drop_na(
         missingness,
         log2_pg_maxlfq
      ) %>%
      dplyr::group_by(
         comparison,
         protein_group
      ) %>%
      mutate(
         n_obs = dplyr::n()
      ) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(
         protein_group,
         comparison,
         missingness,
         .data$n_obs
      )

   moderated_t_test_result <- purrr::map_dfr(
      .x = names,
      .f = ~ limma::topTable(
         moderated_t_test_fit3,
         coef    = .x,
         number  = Inf,
         confint = T,
         sort.by = 'p'
      ) %>%
         tibble::rownames_to_column(
            'protein_group'
         ) %>%
         dplyr::mutate(
            comparison = .x
         )
   ) %>%
      dplyr::mutate(
         comparison     = stringr::str_replace_all(
            comparison,
            pattern     = '^x|(?<=_vs_)x',
            replacement = ' '
         )
      ) %>%
      dplyr::rename(
         diff          = .data$logFC,
         CI_2.5        = .data$CI.L,
         CI_97.5       = .data$CI.R,
         t_statistic   = .data$t,
         avg_abundance = .data$AveExpr,
         pval        = .data$P.Value
      ) %>%
      dplyr::left_join(
         moderated_test_missingness,
         by = c(
            'protein_group',
            'comparison'
         )
      ) %>%
      select(-c(
         adj.P.Val,
         t_statistic,
         B,
         missingness,
         n_obs
      ))

   return(moderated_t_test_result)
   message('DONE', appendLF = F)

}
