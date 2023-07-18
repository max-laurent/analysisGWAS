#' Manhattan Plot with specific markers with zoom
#'
#' Creates a basic Manhattan plot for genome-wide association studies (GWAS) with SNPs in specific colors.
#'
#' @param data Input data containing information for the plot. From GNPasso.
#' @param pos_markers Name of the column in the `data` that contains marker positions.
#' @param name_markers The name of the column in the `data` that contains marker names.
#' @param y Name of the column in the `data` that contains the y-axis values (usually -log10(p-value)).
#' @param chr_col Name of the column in the `data` that contains chromosome information.
#' @param color_col Name of the column containing the different categories to color the SNPs.
#' @param save_time A logical value indicating whether to discard 95% of the SNP that are not significant
#' @param GWAS_name The name of the GWAS being conducted.
#' @param discard_chr The chromosome to be discarded from the plot.
#' @param list_keep Vector containing the name of the markers to keep in the filtering process even if they are not significant
#' @param MAF_filter A minor allele frequency filter value.
#' @param MAF_col Name of the column in the `data` that contains minor allele frequency information.
#' @param chr Int of the chromosome where the markers you want to plot are
#' @param start Position of the first marker you want to plot
#' @param stop Position of the last marker you want to plot
#'
#' @return Manhattan plot
#' @export
#' @author Maxime Laurent
#' @examples
#' # Example usage:
#' manhattan_plot_color(data, pos_markers = "Marker_Position", y = "LogPval",
#'                     chr_col = "Chromosome", color_chr = c("gray", "darkgray"),
#'                     save_time = TRUE, GWAS_name = NULL, discard_chr = 11,
#'                     MAF_filter = 0.05, MAF_col = "MAF")
manhattan_plot_basic_zoom <- function(data, pos_markers = "Marker_Position", name_markers = "Marker_Name", y = "LogPval", chr_col = "Chromosome", GWAS_name = NULL, discard_chr = 11, MAF_filter = 0.05, MAF_col = "MAF", save_time = TRUE, list_keep = NULL, chr = NULL, start = NULL, stop = NULL){

  data <- process_run(data = data, pos_markers = pos_markers, y = y, chr_col = chr_col, save_time = save_time, discard_chr = discard_chr, MAF_filter = MAF_filter, MAF_col = MAF_col, GWAS_name = GWAS_name, list_keep = list_keep)
  max_y <- base::ceiling((max(data[[y]]) + 4)/ 5 * 5)

  if(is.null(chr)){
    chr <- 1
  }
  if(is.null(start)){
    start <- min(data %>% filter(.data[[chr_col]] == chr) %>% pull(.data[[pos_markers]]))
  }
  if(is.null(stop)){
    stop <- max(data %>% filter(.data[[chr_col]] == chr) %>% pull(.data[[pos_markers]]))
  }


  data <- data %>% filter(.data[[chr_col]] == chr & .data[[pos_markers]] >= start & .data[[pos_markers]] <= stop)


  # Create a list containing the information for the mapping of the data
  mapping <- aes(x = Marker_Position, y= .data[[y]], color = as.factor(.data[[chr_col]]))

  ggplot2::ggplot(data, mapping = mapping) +
    ggplot2::geom_point(alpha = 0.8, size = 1) +
    ggplot2::scale_color_manual(values = rep(c("gray", "darkgray"), length(levels(as.factor(data[[chr_col]]))))) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
    # custom X axis:
    ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
    ggplot2::scale_y_continuous(limits= c(0, max_y), expand = c(0, 0) ) +
    # remove space between plot area and x axis
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](pvalue))) +
    # Custom the theme:
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line.y = element_line(colour = "black"),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid = element_blank(),
      #panel.grid.minor.x = element_blank()
    )

}








