#' scale to (0,1) interval
#' @description scale to (0,1) interval
#' @export
#' @param x numeric vector
#' @return This function returns the vector \code{x} scaled to lie within the (0,1) interval


scale_01 <- function(x){
  (x - min(x))/max(x - min(x))
}
