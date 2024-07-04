#' @title get_smoothPara
#' @description
#' Get smooth parameters list.
#'
#' @param smooth Whether to perform smooth.
#' @param method Smooth method, sg or mean.
#' @param size SmoothMean size.
#' @param p sg p.
#' @param n sg n.
#' @param m sg m.
#' @param ts sg ts.
#'
#' @return A parameter list.
#' @export
#'
#' @examples
#' smoothPara <- get_smoothPara()
get_smoothPara <- function(smooth = TRUE, method = "mean", size = 3,
                       p = 3, n = p + 3 - p%%2, m = 0, ts = 1){
  if(method != "mean" & method != "sg") stop("method must be mean or sg!")
  if(size %% 2 == 0){
    stop("size should be a singular!")
  }
  return(list(smooth = smooth, method = method,
              size = size,
              p = p, n = n, m = m, ts = 1))
}
#' @title get_baselinePara
#' @description
#' Get baselineEs function parameters.
#'
#' @param threshold threshold.
#' @param tol_m tol_m.
#' @param loops loops.
#'
#' @return A parameters list.
#' @export
#'
#' @examples
#' baselinePara <- getbaselinePara()
get_baselinePara <- function(threshold = 1, tol_m = 30, loops = 6){
  return(list(threshold = threshold, tol_m = tol_m, loops = loops))
}
