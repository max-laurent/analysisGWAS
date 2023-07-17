#' Process for manhattan plot
#'
#' Processes the input data for further analysis in a GWAS run. To prepare the data for a Manhattan plot
#'
#' @param data output from GNPasso
#' @param pos_markers The name of the column in the `data` that contains marker positions.
#' @param name_markers The name of the column in the `data` that contains marker names.
#' @param y The name of the column in the `data` that contains the y-axis values (usually -log10(p-value)).
#' @param chr_col The name of the column in the `data` that contains chromosome information.
#' @param GWAS_name The name of the GWAS being conducted.
#' @param save_time A logical value indicating whether to discard 95% of the SNP that are not significant
#' @param list_keep Vecotr containing the name of the markers to keep in the filtering process even if they are not significant
#' @param discard_chr The chromosome to be discarded from the plot.
#' @param MAF_filter A minor allele frequency filter value.
#' @param MAF_col Name of the column in the `data` that contains minor allele frequency information.
#'
#' @return processed data
#' @export
#'
#' @author Maxime Laurent
#'
#' @examples
#' # Example usage:
#' processed_data <- process_run(data, pos_markers = "Marker_Position", name_markers = "Marker_Name",
#'                               y = "LogPval", chr_col = "Chromosome", GWAS_name = NULL,
#'                               save_time = TRUE, discard_chr = 11, MAF_filter = 0.05, MAF_col = "MAF")
process_run <- function(data, pos_markers = "Marker_Position", name_markers = "Marker_Name", y = "LogPval", chr_col = "Chromosome", GWAS_name = NULL, save_time = TRUE, list_keep = NULL, discard_chr = 11, MAF_filter = 0.05, MAF_col = "MAF"){
  # Control the function parameters
  assertthat::assert_that(is.numeric(pos_markers)|is.character(pos_markers), msg="pos_markers must be text or int")
  assertthat::assert_that(is.numeric(name_markers)|is.character(name_markers), msg="name_markers must be text or int")
  assertthat::assert_that(is.numeric(y)|is.character(y), msg="y must be text or int")
  assertthat::assert_that(is.numeric(chr_col)|is.character(chr_col), msg="chr_col must be text or int")
  assertthat::assert_that(is.null(GWAS_name)|is.numeric(GWAS_name)|is.character(GWAS_name), msg="GWAS_name must be text or int")
  assertthat::assert_that(is.logical(save_time), msg="save_time must be boolean")
  assertthat::assert_that(is.numeric(discard_chr), msg="discard_chr must be int")
  assertthat::assert_that(is.null(MAF_col)| (MAF_col %in% names(data)), msg="MAF column must be in data columns or NULL")
  assertthat::assert_that(is.null(list_keep)|is.vector(list_keep), msg = "list_keep must be a vector or NULL")

  # Calculate cumulative position based on chromosome and marker positions
  processed_data <- data %>%
    dplyr::group_by(.data[[chr_col]]) %>%
    dplyr::summarise(chr_length = max(.data[[pos_markers]])) %>%
    dplyr::mutate(tot = cumsum(chr_length) - chr_length) %>%
    dplyr::select(-chr_length) %>%
    dplyr::left_join(data,. , by = (chr_col = chr_col)) %>%
    dplyr::arrange(.data[[chr_col]], .data[[pos_markers]]) %>%
    dplyr::mutate(Marker_Position_cum = .data[[pos_markers]] + tot) %>%
    {if(MAF_col %in% names(data)) dplyr::filter(., MAF > 0.05 ) else .} %>%
    {if(!is.null(GWAS_name)) dplyr::mutate(., GWAS = GWAS_name) else .} %>%
    dplyr::filter(.data[[chr_col]] != 11)

  # Filter 95% of the non significant markers
  if(save_time){
    th <- base::sample(processed_data[processed_data[[y]] < 5, name_markers],
                 round(0.05 * length(processed_data[processed_data[[y]] < 5, name_markers])))
    processed_data <- processed_data %>% dplyr::filter(.data[[y]] >= 5 | .data[[name_markers]] %in% th | .data[[name_markers]] %in% list_keep)
  }
  return(processed_data)
}

#' Calculate the axis of the manhattan plot
#'
#' @param data the processed data
#' @param chr_col column containing the chromosome position
#'
#' @return dataframe containing the axis
#' @export
#' @author Maxime Laurent
#' @examples
calc_axis <- function(data, chr_col = "Chromosome"){
  axis <- data %>%
    dplyr::group_by(.data[[chr_col]]) %>%
    dplyr::summarise(max = max(Marker_Position_cum), min = min(Marker_Position_cum)) %>%
    dplyr::mutate(center = (max/2 + min/2))
  return(axis)
}


