#' Render isobands
#'
#' This function generates a grid grob that represents isobands.
#'
#' @param bands Isobands, as produced by the [`isobands()`] function.
#' @param gp Grid graphical parameters. Parameters are recycled among
#'   the total number of bands drawn.
#' @param units A character string specifying the units in which to
#'   interpret the isobands coordinates. Defaults to `"npc"`.
#' @seealso
#' See [`isolines_grob()`] for drawing of isolines.
#' @examples
#' library(grid)
#'
#' viridis_pal <- colorRampPalette(
#'   c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725"),
#'   space = "Lab"
#' )
#'
#' x <- (1:ncol(volcano))/(ncol(volcano)+1)
#' y <- (nrow(volcano):1)/(nrow(volcano)+1)
#' bands <- isobands(x, y, volcano, 5*(18:38), 5*(19:39))
#'
#' b <- isobands_grob(
#'   bands,
#'   gp = gpar(col = "black", fill = viridis_pal(21), alpha = 0.5)
#' )
#'
#' grid.newpage()
#' grid.draw(b)
#' @export
isobands_grob <- function(bands, gp = gpar(), units = "npc") {
  gTree(
    bands = bands,
    gp_user = gp,
    units = units,
    cl = "isobands_grob"
  )
}

#' @export
makeContent.isobands_grob <- function(x) {
  make_bands_grobs <- function(data, col, fill, alpha, lty, lwd, lex, lineend, linejoin, linemitre) {
    if (length(data$x) == 0) {
      return(NULL)
    }
    pathGrob(
      data$x, data$y, data$id,
      default.units = x$units,
      gp = gpar(
        col = col, fill = fill, alpha = alpha, lty = lty, lwd = lwd, lex = lex,
        lineend = lineend, linejoin = linejoin, linemitre = linemitre
      )
    )
  }

  # merge current and grob-specific graphical parameters so we can redistribute among isolevels
  gp <- modifyList(get.gpar(), x$gp_user)
  n <- length(x$bands)

  bands_grobs <- mapply(
    make_bands_grobs,
    x$bands,
    rep_len(gp$col, n),
    rep_len(gp$fill, n),
    rep_len(gp$alpha, n),
    rep_len(gp$lty, n),
    rep_len(gp$lwd, n),
    rep_len(gp$lex, n),
    rep_len(gp$lineend, n),
    rep_len(gp$linejoin, n),
    rep_len(gp$linemitre, n),
    SIMPLIFY = FALSE
  )

  children <- do.call(gList, bands_grobs)
  setChildren(x, children)
}
