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
#'
#' m <- volcano
#' b <- isobands((1:ncol(m))/(ncol(m)+1), (nrow(m):1)/(nrow(m)+1), m,
#'               10*9:19, 10*10:20)
#' bands <- iso_to_sfg(b)
#' x <- st_sf(level = as.numeric(sub(":.*", "", names(bands))), geometry = st_sfc(bands))
#' ggplot(x) + geom_sf(aes(color = level, fill = level))
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
  mapply(multilinestring, x, SIMPLIFY = FALSE)
}

multilinestring <- function(object) {
  x <- split(object$x, object$id)
  y <- split(object$y, object$id)
  structure(
    unname(mapply(cbind, x, y, SIMPLIFY = FALSE)),
    class = c("XY", "MULTILINESTRING", "sfg")
  )
}

#' @export
iso_to_sfg.isobands <- function(x) {
  mapply(multipolygon, x, SIMPLIFY = FALSE)
}

multipolygon <- function(object) {
  separate_polygons(object$x, object$y, object$id)
}

