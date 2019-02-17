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
#' @param units A character string specifying the units in which to
#'   interpret the isolines coordinates. Defaults to `"npc"`.
#' @param label_col Color applied to labels. Can be used to override the
#'   color provided in `gp`, in case labels and lines should have different
#'   colors.
#' @param label_alpha Alpha applied to labels. Can be used to override the
#'   alpha value provided in `gp`, in case labels and lines should have
#'   different alpha values.
#' @examples
#' x <- (1:ncol(volcano))/(ncol(volcano)+1)
#' y <- (nrow(volcano):1)/(nrow(volcano)+1)
#' lines <- isolines(x, y, volcano, 10*(10:18))
#'
#' library(grid)
#' grid.newpage()
#'
#' # create a colored background
#' grid.newpage()
#' grid.rect(width = .8, height = .8, gp = gpar(fill = "azure1", col = NA))
#' grid.rect(width = .6, height = .6, gp = gpar(fill = "cornsilk1", col = NA))
#' grid.rect(width = .4, height = .4, gp = gpar(fill = "cornsilk2", col = NA))
#' grid.rect(width = .2, height = .2, gp = gpar(fill = "cornsilk3", col = NA))
#'
#' # draw labeled lines
#' grid.draw(isolines_grob(lines))
#' @export
isolines_grob <- function(lines, gp = gpar(), breaks = NULL, labels = NULL,
                          units = "npc", label_col = NULL, label_alpha = NULL) {
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

  rows <- mapply(
    place_labels,
    lines[breaks],
    breaks,
    labels,
    match(breaks, names(lines)), # index of labeled lines in original list of lines, for matching of graphical parameters
    1:length(breaks), # index into original list of breaks, for matching graphical parameters
    SIMPLIFY = FALSE
  )
  labels_data <- Reduce(rbind, rows)

  gTree(
    lines = lines,
    breaks = breaks,
    labels = labels,
    labels_data = labels_data,
    gp_combined = gp,
    label_col = label_col,
    label_alpha = label_alpha,
    units = units,
    cl = "isolines_grob"
  )
}

#' @export
makeContext.isolines_grob <- function(x) {
  # store current graphics parameters for later
  x$gp_cur <- get.gpar()

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

  # calculate label widths and heights in npc units
  label_widths <- convertWidth(stringWidth(labels_data$label), x$units, valueOnly = TRUE)
  label_heights <- convertHeight(
    stringHeight(labels_data$label) + stringDescent(labels_data$label),
    x$units, valueOnly = TRUE
  )

  # calculate the clip box for each label
  clip_boxes <- data.frame(
    x = labels_data$x, y = labels_data$y,
    width = label_widths, height = label_heights,
    theta = labels_data$theta
  )

  make_lines_grobs <- function(data, col, alpha, lty, lwd, lex, lineend, linejoin, linemitre) {
    if (length(data$x) == 0) {
      return(NULL)
    }
    clipped <- clip_lines(data$x, data$y, data$id, clip_boxes)
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

  # calculate rotation angles for text labels
  rot <- 360*labels_data$theta/(2*pi)
  rot <- ifelse(rot <= -90, 180 + rot, rot)
  rot <- ifelse(rot > 90, rot - 180, rot)

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
    labels_data$label, labels_data$x, labels_data$y, rot = rot,
    default.units = x$units,
    gp = gpar(col = col, alpha = alpha)
  )

  children <- do.call(gList, c(lines_grobs, list(labels_grob)))
  setChildren(x, children)
}


place_labels <- function(line_data, break_id, label, index, break_index) {
  # return empty row if either missing line data or missing label
  if (length(line_data$x) == 0 || is.na(label)) {
    return(
      data.frame(
        index = integer(0),
        break_index = integer(0),
        break_id = character(0), label = character(0),
        x = numeric(0), y = numeric(0), theta = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }

  # find location for labels
  idx <- c(
    which(line_data$y == min(line_data$y))[1],
    which(line_data$y == max(line_data$y))[1]
  )

  # calculate label position
  pos <- label_positions(line_data, idx)

  # return results
  data.frame(
    index = index, # index into original list of all isolines
    break_index = break_index, # index into original list of breaks
    break_id = break_id,
    label = label,
    x = pos$x, y = pos$y, theta = pos$theta,
    stringsAsFactors = FALSE
  )
}

# Calculate the position and rotation of labels based
# on the x, y data and the index position(s) in the data.
# The variable `n` sets the neighborhood size, n = 2 means
# two points in either direction are used.
label_positions <- function(data, idx, n = 2) {
  out <- data.frame(x = numeric(0), y = numeric(0), theta = numeric(0))

  for (i in seq_along(idx)) {
    out[i, ] <- label_position_single(data, idx[i], n)
  }

  out
}

# function that does all the work for label_positions().
# requires a single index idx
label_position_single <- function(data, idx, n = 2) {
  # set of indices belonging to this label
  idx_set <- which(data$id == data$id[idx])
  idx_min <- min(idx_set)
  idx_max <- max(idx_set)

  # if the first and the last point are the same we wrap, otherwise we truncate
  if (data$x[idx_min] == data$x[idx_max] && data$y[idx_min] == data$y[idx_max]) {
    idx_range <- (idx_max - idx_min)
    i <- ((idx - n):(idx + n)-idx_min) %% idx_range + idx_min
  } else {
    i <- (max(idx - n, idx_min):min(idx + n, idx_max))
  }

  x <- data$x[i]
  y <- data$y[i]
  xave <- mean(x)
  yave <- mean(y)
  m <- cbind(x - xave, y - yave)
  v <- svd(m)$v
  list(x = xave, y = yave, theta = atan2(v[2], v[1]))
}
