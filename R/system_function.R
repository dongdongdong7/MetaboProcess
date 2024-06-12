#' @title get_os
#' @description
#' Get operating system type.
#'
#' @return A character.
#' @export
#'
#' @examples
#' get_os()
get_os <- function(){
  return(base::Sys.info()[["sysname"]])
}
