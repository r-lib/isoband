#' Draw polygon path legend keys
#'
#' @usage NULL
#' @format NULL
#' @keywords internal
#' @export
draw_key_polypath <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)

  idx <- is.null(data$fill_alpha)
  data$fill_alpha[idx] <- data$alpha[idx]

  grobTree(
    rectGrob(
      width = unit(1, "npc") - unit(lwd, "mm"),
      height = unit(1, "npc") - unit(lwd, "mm"),
      gp = gpar(
        col = NA,
        fill = alpha(data$fill, data$fill_alpha)
      )
    ),
    rectGrob(
      width = unit(1, "npc") - unit(lwd, "mm"),
      height = unit(1, "npc") - unit(lwd, "mm"),
      gp = gpar(
        col = alpha(data$colour, data$alpha),
        fill = NA,
        lty = data$linetype,
        lwd = lwd * .pt,
        linejoin = "mitre"
      )
    )
  )
}
