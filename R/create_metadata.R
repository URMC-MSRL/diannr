#' Create Metadata Information
#'
#' @description Creates metadata information, associating treatment groups to sample names
#' @param data The data frame output from the prepare_data function
#'
#' @return A matrix that has a condition name for every sample name
#' @export
#'
create_metadata <- function(data) {

   #Pull out just the sample names from data
   sample_names <- data %>%
      dplyr::distinct(sample)
   conditions <- readline('What is the condition of each sample in order? (ensure to seaprate each condition by a comma)')
   conditions <- strsplit(conditions, ',')

   sample_annotation <- data.frame(
      sample_names,
      conditions
   ) %>%
      purrr::set_names(
         'sample',
         'condition'
      )
}
