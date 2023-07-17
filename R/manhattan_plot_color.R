#' Manhattan Plot with specific markers in color
#'
#' Creates a basic Manhattan plot for genome-wide association studies (GWAS).
#'
#' @param data Input data containing information for the plot. From GNPasso.
#' @param pos_markers Name of the column in the `data` that contains marker positions.
#' @param name_markers The name of the column in the `data` that contains marker names.
#' @param y Name of the column in the `data` that contains the y-axis values (usually -log10(p-value)).
#' @param chr_col Name of the column in the `data` that contains chromosome information.
#' @param color_col A vector specifying the colors to be used for different SNPs.
#' @param save_time A logical value indicating whether to discard 95% of the SNP that are not significant
#' @param GWAS_name The name of the GWAS being conducted.
#' @param discard_chr The chromosome to be discarded from the plot.
#' @param list_keep Vector containing the name of the markers to keep in the filtering process even if they are not significant
#' @param MAF_filter A minor allele frequency filter value.
#' @param MAF_col Name of the column in the `data` that contains minor allele frequency information.
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
manhattan_plot_color <- function(data, pos_markers = "Marker_Position", name_markers = "Marker_Name", y = "LogPval", chr_col = "Chromosome", color_col =  NULL, save_time = TRUE, GWAS_name = NULL, discard_chr = 11, MAF_filter = 0.05, MAF_col = "MAF", list_keep = NULL){

  data <- process_run(data = data, pos_markers = pos_markers, y = y, chr_col = chr_col, save_time = save_time, discard_chr = discard_chr, MAF_filter = MAF_filter, MAF_col = MAF_col, GWAS_name = GWAS_name, list_keep = list_keep)
  axis <- calc_axis(data = data, chr_col = "Chromosome")
  ncol <- length(unique(data[[color_col]]))-2
  data[[color_col]] <- as.factor(data[[color_col]])
  max_y <- base::ceiling((max(data[[y]]) + 4)/ 5 * 5)
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")

  ggplot2::ggplot(data, aes(x = Marker_Position_cum, y= .data[[y]])) +
    ggplot2::geom_point(aes(color = .data[[color_col]]), alpha = 0.8, size = 0.4) +
    ggplot2::scale_color_manual(values = c(pal[1:ncol], "gray", "darkgray")) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
    # custom X axis:
    ggplot2::scale_x_continuous(label = axis[[chr_col]], breaks= axis$center, expand = c(0.01,0.01)) +
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








