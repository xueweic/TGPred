
#' Calculate the gridient of Huber function for an input value z.
#'
#' @param z Input value z.
#' @param M Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).
#'
#' @return The gridient of Huber function for an input value z
#' @export
#'
#' @examples
Huber_grid <- function(z, M) {
  ifelse(abs(z) < M, 2 * z, 2 * M * sign(z))
}

#' Calculate the Huber function for an input value z.
#'
#' @param z Input value z
#' @param M Shape parameter, which is defaulted to be one-tenth of the interquartile range (IRQ).
#'
#' @return The Huber function for an input value z
#' @export
#'
#' @examples
HuberFunction <- function(z, M) {
  ifelse(abs(z) < M, z^2, 2 * M * abs(z) - M^2)
}
