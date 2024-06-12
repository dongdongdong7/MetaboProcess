#' @title plot_chrDf
#' @description
#' Plot a chromatogram using chromatographic data frame.
#'
#' @param chrDf chrDf.
#' @param linewidth linewidth.
#' @param noise noise.
#' @param xlim xlim.
#' @param baseline
#' Whether to draw the baseline, if TRUE, you need to ensure that there is a baseline in chrDf.
#'
#' @return ggplot object.
#' @export
#'
#' @examples
#' plot_chrDf(chrDf = chrDf_i, linewidth = 2)
plot_chrDf <- function(chrDf, linewidth = 1, noise = NA, xlim = NA, baseline = FALSE){
  if(is.numeric(xlim)){
    chrDf <- chrDf %>%
      dplyr::filter(rt >= xlim[1] & rt <= xlim[2])
    ymax <- max(chrDf$intensity)
  }
  title <- paste0(round(min(chrDf$mz), 4), " - ", round(max(chrDf$mz), 4))
  p <- ggplot2::ggplot(data = chrDf, mapping = ggplot2::aes(x = rt)) +
    ggplot2::geom_line(mapping = ggplot2::aes(y = intensity), linewidth = linewidth) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Retention Time", y = "Intensity", title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 15),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text = ggplot2::element_text(size = 10))
  if(!is.na(noise)){
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept = noise),
                                 color = "#990000",
                                 linewidth = linewidth,
                                 linetype = "dashed")
  }
  if(is.numeric(xlim)){
    p <- p + ggplot2::ylim(c(0, ymax))
  }
  if(baseline){
    p <- p + ggplot2::geom_line(ggplot2::aes(x = rt, y = baseline), col = "red", linewidth = 1, linetype = "dashed")
  }
  return(p)
}
