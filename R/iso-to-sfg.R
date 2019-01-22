#' Convert isolines or isobands to sfg object
#'
#' Convert isolines or isobands to an sf geometry collection (`sfg`) object. Further downstream
#' processing needs to happen via the sf package.
#' @param x The object to convert.
#' @examples
#' library(sf)
#' library(ggplot2)
#'
#' m <- matrix(c(0, 2, 2, 2, 0,
#'               0, 1, 0, 1, 0,
#'               0, 1, 0, 0, 0,
#'               0, 1, 0, 1, 0,
#'               0, 0, 0, 0, 0), 5, 5, byrow = TRUE)
#'
#' z <- isolines(1:ncol(m), nrow(m):1, m, c(0.5, 1.5))
#' lines <- iso_to_sfg(z)
#' x <- st_sf(level = names(lines), geometry = st_sfc(lines))
#' ggplot(x) + geom_sf(aes(color = level))
#' @export
iso_to_sfg <- function(x) {
  UseMethod("iso_to_sfg", x)
}

#' @export
iso_to_sfg.default <- function(x) {
  stop(
    "Cannot convert objects of type ", paste(class(x), collapse = "/"), " to sf.",
    call. = FALSE
  )
}

#' @export
iso_to_sfg.isolines <- function(x) {
  mapply(multilinestring, x)
}

multilinestring <- function(object) {
  x <- split(object$x, object$id)
  y <- split(object$y, object$id)
  structure(
    unname(mapply(cbind, x, y, SIMPLIFY = FALSE)),
    class = c("XY", "MULTILINESTRING", "sfg")
  )
}
