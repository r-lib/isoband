#' Render labeled isolines
#'
#' This function generates a grid grob that represents labeled isolines.
#'
#' @param lines Isolines, as produced by the [`isolines()`] function.
#' @param gp Grid graphical parameters. Parameters applying to lines
#'   (such as `col`, `lwd`, `lty`, etc.) are recycled among the total
#'   number of lines drawn. Parameters applying only to labels (such
#'   as `fontfamily`, `fontsize`) are recycled among the specified
#'   breaks only. The two parameters `col` and `alpha` are also applied
#'   to labels, unless overridden (see `label_col` and `label_alpha`),
#'   but are matched to the corresponding lines.
#' @param breaks Character vector specifying the isolines that should be
#'   labeled. If `NULL`, labels all isolines.
#' @param labels Character vector specifying the labels for each break.
#'   If `NULL`, uses the breaks as labels. The number of labels provided
#'   must match the number of breaks provided.
#' @param margin Unit object of length 4 specifying the top, right, bottom,
#'   and left margins around each text label. The same margins are applied
#'   to all labels.
#' @param label_col Color applied to labels. Can be used to override the
#'   color provided in `gp`, in case labels and lines should have different
#'   colors.
#' @param label_alpha Alpha applied to labels. Can be used to override the
#'   alpha value provided in `gp`, in case labels and lines should have
#'   different alpha values.
#' @param label_placer Function that controls how labels are placed along
#'   the isolines. Uses [`label_placer_minmax()`] by default.
#' @param units A character string specifying the units in which to
#'   interpret the isolines coordinates. Defaults to `"npc"`.
#' @seealso
#' See [`isobands_grob()`] for drawing of isobands. See [`label_placer_minmax()`] for
#' label placement strategies.
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
#' lines <- isolines(x, y, volcano, 5*(19:38))
#' bands <- isobands(x, y, volcano, 5*(18:38), 5*(19:39))
#'
#' b <- isobands_grob(
#'   bands,
#'   gp = gpar(col = NA, fill = viridis_pal(21), alpha = 0.4)
#' )
#' l <- isolines_grob(
#'   lines, breaks = 20*(5:10),
#'   gp = gpar(
#'     lwd = c(.3, 1, .3, .3)
#'   )
#' )
#'
#' grid.newpage()
#' grid.draw(b)
#' grid.draw(l)
#' @export
isolines_grob <- function(lines, gp = gpar(), breaks = NULL, labels = NULL,
                          margin = unit(c(1, 1, 1, 1), "pt"), label_col = NULL, label_alpha = NULL,
                          label_placer = label_placer_minmax(), units = "npc") {
  if (is.null(breaks)) {
    breaks <- names(lines)
  } else {
    breaks <- as.character(breaks)
  }

  if (is.null(labels)) {
    labels <- breaks
  } else if (length(labels) != length(breaks)) {
    stop("Number of labels must match the number of breaks.", call. = FALSE)
  } else {
    labels <- as.character(labels)
  }

  if (length(margin) != 4 || !is.unit(margin)) {
    stop("The `margin` parameter must be a unit object of length four.", call. = FALSE)
  }

  # first set up a data frame with all the label info
  labels_data <- data.frame(
    index = match(breaks, names(lines)), # index of labeled lines in original list of lines, for matching of graphical parameters
    break_index = 1:length(breaks), # index into original list of breaks, for matching graphical parameters
    break_id = breaks, # identifier for each break (corresponds to the name in the lines column)
    label = labels, # label for each break
    stringsAsFactors = FALSE
  )
  # then calculate the position of all labels via the `label_placer()` function
  labels_data <- label_placer(lines, labels_data)

  gTree(
    lines = lines,
    breaks = breaks,
    labels = labels,
    labels_data = labels_data,
    gp_combined = gp,
    margin = margin,
    label_col = label_col,
    label_alpha = label_alpha,
    label_placer = label_placer,
    units = units,
    cl = "isolines_grob"
  )
}

#' @export
makeContext.isolines_grob <- function(x) {
  # store current graphics parameters for later
  x$gp_cur <- get.gpar()

  # if we have no labels then nothing else needs to be done
  if (nrow(x$labels_data) == 0) {
    return(x)
  }

  # we need to set up the gp slot for text labels, so font calculations work
  gp <- modifyList(x$gp_cur, x$gp_combined)

  # map graphical parameters to duplicated labels
  # we only handle font parameters here, to guarantee proper
  # calculation of label sizes
  n <- length(x$breaks)
  cex <- rep_len(gp$cex, n)[x$labels_data$break_index]
  lineheight <- rep_len(gp$lineheight, n)[x$labels_data$break_index]
  fontfamily <- rep_len(gp$fontfamily, n)[x$labels_data$break_index]
  fontsize <- rep_len(gp$fontsize, n)[x$labels_data$break_index]

  # fontface needs special treatment, since it can be NULL if font
  # is specified
  if (is.null(gp$fontface)) {
    # we work with font
    font <- rep_len(gp$font, n)[x$labels_data$break_index]

    x$gp <- gpar(
      cex = cex, fontsize = fontsize, lineheight = lineheight,
      fontfamily = fontfamily, font = font
    )
  } else {
    # we work with fontface
    fontface <- rep_len(gp$fontface, n)[x$labels_data$break_index]

    x$gp <- gpar(
      cex = cex, fontsize = fontsize, lineheight = lineheight,
      fontfamily = fontfamily, fontface = fontface
    )
  }

  x
}

#' @export
makeContent.isolines_grob <- function(x) {
  labels_data <- x$labels_data

  if (nrow(labels_data) == 0) {
    # no labels to be drawn, nothing to be done
    return(isolines_grob_makeContent_nolabels(x))
  }

  # calculate label widths and heights in npc units
  label_widths <- convertWidth(stringWidth(labels_data$label), x$units, valueOnly = TRUE)
  label_heights <- convertHeight(
    stringHeight(labels_data$label) + stringDescent(labels_data$label),
    x$units, valueOnly = TRUE
  )

  # get viewport aspect ratio to correct clipping for rotated labels
  asp <- convertHeight(unit(1, "pt"), x$units, valueOnly = TRUE) / convertWidth(unit(1, "pt"), x$units, valueOnly = TRUE)
  #print(asp)

  # calculate margins in npc units
  margin_rl <- convertWidth(x$margin[c(2, 4)], x$units, valueOnly = TRUE)
  margin_tb <- convertHeight(x$margin[c(1, 3)], x$units, valueOnly = TRUE)
  margin_w <- sum(margin_rl)
  margin_h <- sum(margin_tb)
  margin_wdiff <- margin_rl[2] - margin_rl[1]
  margin_hdiff <- margin_rl[1] - margin_rl[2]

  # calculate the clip box for each label
  # xoff and yoff are needed to correct for uneven label margins
  xoff <- -margin_wdiff*cos(labels_data$theta)/2 + (margin_hdiff/asp)*sin(labels_data$theta)/2
  yoff <- -margin_wdiff*sin(labels_data$theta)/2 - margin_hdiff*cos(labels_data$theta)/2

  clip_boxes <- data.frame(
    x = labels_data$x + xoff, y = labels_data$y + yoff,
    width = label_widths + margin_w, height = label_heights + margin_h,
    theta = labels_data$theta
  )

  make_lines_grobs <- function(data, col, alpha, lty, lwd, lex, lineend, linejoin, linemitre) {
    if (length(data$x) == 0) {
      return(NULL)
    }
    clipped <- clip_lines(data$x, data$y, data$id, clip_boxes, asp = asp)
    if (length(clipped$x) == 0) {
      return(NULL)
    }
    polylineGrob(
      clipped$x, clipped$y, clipped$id,
      default.units = x$units,
      gp = gpar(
        col = col, alpha = alpha, lty = lty, lwd = lwd, lex = lex,
        lineend = lineend, linejoin = linejoin, linemitre = linemitre
      )
    )
  }

  # merge current and grob-specific graphical parameters so we can redistribute among isolevels
  gp <- modifyList(x$gp_cur, x$gp_combined)
  n <- length(x$lines)

  lines_grobs <- mapply(
    make_lines_grobs,
    x$lines,
    rep_len(gp$col, n),
    rep_len(gp$alpha, n),
    rep_len(gp$lty, n),
    rep_len(gp$lwd, n),
    rep_len(gp$lex, n),
    rep_len(gp$lineend, n),
    rep_len(gp$linejoin, n),
    rep_len(gp$linemitre, n),
    SIMPLIFY = FALSE
  )

  # calculate color and alpha for text labels
  if (is.null(x$label_col)) {
    col <- rep_len(gp$col, length(x$lines))[x$labels_data$index]
  } else {
    col <- rep_len(x$label_col, length(x$breaks))[x$labels_data$break_index]
  }
  if (is.null(x$label_alpha)) {
    # alpha is cumulative, so we don't take it from the current viewport
    alpha <- rep_len(x$gp_combined$alpha %||% 1, length(x$lines))[x$labels_data$index]
  } else {
    alpha <- rep_len(x$label_alpha, n)[x$labels_data$break_index]
  }

  labels_grob <- textGrob(
    labels_data$label, labels_data$x, labels_data$y, rot = 360*labels_data$theta/(2*pi),
    default.units = x$units,
    gp = gpar(col = col, alpha = alpha)
  )

  children <- do.call(gList, c(lines_grobs, list(labels_grob)))
  setChildren(x, children)
}

isolines_grob_makeContent_nolabels <- function(x) {
  make_lines_grobs <- function(data, col, alpha, lty, lwd, lex, lineend, linejoin, linemitre) {
    if (length(data$x) == 0) {
      return(NULL)
    }

    polylineGrob(
      data$x, data$y, data$id,
      default.units = x$units,
      gp = gpar(
        col = col, alpha = alpha, lty = lty, lwd = lwd, lex = lex,
        lineend = lineend, linejoin = linejoin, linemitre = linemitre
      )
    )
  }

  # merge current and grob-specific graphical parameters so we can redistribute among isolevels
  gp <- modifyList(x$gp_cur, x$gp_combined)
  n <- length(x$lines)

  lines_grobs <- mapply(
    make_lines_grobs,
    x$lines,
    rep_len(gp$col, n),
    rep_len(gp$alpha, n),
    rep_len(gp$lty, n),
    rep_len(gp$lwd, n),
    rep_len(gp$lex, n),
    rep_len(gp$lineend, n),
    rep_len(gp$linejoin, n),
    rep_len(gp$linemitre, n),
    SIMPLIFY = FALSE
  )

  children <- do.call(gList, lines_grobs)
  setChildren(x, children)
}

