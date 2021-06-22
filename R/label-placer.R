#' Generic label placement function
#'
#' The simple label placer processes separate isolines independently and places
#' labels for each line using a placer function that does the actual placement work.
#' This label placer is not meant to be used by end users, but rather facilitates the
#' development of new label placers, such as [`label_placer_minmax()`].
#' @param lines Isolines object for which labels should be placed.
#' @param labels_data A data frame containing information about which labels should
#'   be placed.
#' @param placer_fun A function that takes an individual isoline plus its associated
#'   break id as input and returns a data frame specifying label positions. The data
#'   frame should have three columns called `x`, `y`, and `theta`. `x` and `y` specify
#'   the label position, and `theta` specifies the label angle in radians. The data
#'   frame can have multiple rows, which results in the same label being placed in
#'   multiple locations.
#' @keywords internal
#' @export
label_placer_simple <- function(lines, labels_data, placer_fun) {
  # Calculate the label position for one set of isolines (one level).
  #
  # The line data is specified as a list of x, y, id. The parameters `index`, `break_index`,
  # `break_id`, and `label` are provided simply so they can be added to the resulting data
  # frame holding label positions
  place_labels_impl <- function(line_data, index, break_index, break_id, label, placer_fun) {
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

    # calculate label position
    pos <- placer_fun(line_data, break_id)

    # return results
    if (nrow(pos) > 0) {
      data.frame(
        index = index,
        break_index = break_index,
        break_id = break_id,
        label = label,
        x = pos$x, y = pos$y, theta = pos$theta,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        index = integer(0),
        break_index = integer(0),
        break_id = character(0),
        label = character(0),
        x = numeric(0), y = numeric(0), theta = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }

  rows <- mapply(
    place_labels_impl,
    lines[labels_data$index],
    labels_data$index, # index of labeled lines in original list of lines, for matching of graphical parameters
    labels_data$break_index, # index into original list of breaks, for matching graphical parameters
    labels_data$break_id,
    labels_data$label,
    MoreArgs = list(placer_fun = placer_fun),
    SIMPLIFY = FALSE
  )
  Reduce(rbind, rows)
}



#' Set up a label placement strategy
#'
#' These functions set up various label placement strategies.
#'
#' `label_placer_minmax()` places labels at the horizontal or vertical minima or maxima of
#' the respective isolines.
#'
#' `label_placer_none()` places no labels at all.
#'
#' `label_placer_manual()` places labels at manually defined locations.
#' 
#' `label_placer_middle()` places labels at the middle of each isoline.
#'
#' @param placement String consisting of any combination of the letters
#'   "t", "r", "b", "l" indicating the placement of labels at the top,
#'   to the right, at the bottom, to the left of the isoline.
#' @param rot_adjuster Function that standardizes the rotation angles of the labels.
#'   See e.g. [`angle_halfcircle_bottom()`].
#' @param n Size of the point neighborhood over which the rotation angle should be
#'   calculated.
#' @rdname label_placer
#' @export
label_placer_minmax <- function(placement = "tb", rot_adjuster = angle_halfcircle_bottom(), n = 2) {
  force_all(placement, rot_adjuster, n)

  placer_fun <- function(line_data, ...) {
    # find location for labels
    idx <- stats::na.omit(
      c(
        which( # placement "top"
          isTRUE(grepl("t", placement, fixed = TRUE)) & line_data$y == max(line_data$y)
        )[1],
        which( # placement "bottom"
          isTRUE(grepl("b", placement, fixed = TRUE)) & line_data$y == min(line_data$y)
        )[1],
        which( # placement "left"
          isTRUE(grepl("l", placement, fixed = TRUE)) & line_data$x == min(line_data$x)
        )[1],
        which( # placement "right"
          isTRUE(grepl("r", placement, fixed = TRUE)) & line_data$x == max(line_data$x)
        )[1]
      )
    )

    out <- data.frame(x = numeric(0), y = numeric(0), theta = numeric(0))

    for (i in seq_along(idx)) {
      out[i, ] <- minmax_impl(line_data, idx[i], n)
    }

    # standardize rotation angles for text labels
    out$theta <- rot_adjuster(out$theta)

    out
  }

  # final placer function
  function(lines, labels_data) {
    label_placer_simple(lines, labels_data, placer_fun)
  }
}

# function that does all the work for the minmax label placer.
# requires a single index idx
minmax_impl <- function(data, idx, n) {
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

#' @rdname label_placer
#' @export
label_placer_none <- function() {
  function(...) {
    data.frame(
      index = integer(0),
      break_index = integer(0),
      break_id = character(0), label = character(0),
      x = numeric(0), y = numeric(0), theta = numeric(0),
      stringsAsFactors = FALSE
    )
  }
}

#' @param breaks Character vector specifying the isolines to be labeled,
#'   as in [`isolines_grob()`].
#' @param x,y,theta Numeric vectors specifying the x and y positions and
#'   angles (in radians) for each label corresponding to each break.
#' @rdname label_placer
#' @export
label_placer_manual <- function(breaks, x, y, theta) {
  # recycle all inputs to the same length
  # also has the side effect of forcing them
  n <- max(length(breaks), length(x), length(y), length(theta))
  breaks <- rep_len(breaks, n)
  x <- rep_len(x, n)
  y <- rep_len(y, n)
  theta <- rep_len(theta, n)

  placer_fun <- function(line_data, break_id) {
    idx <- (breaks == break_id)
    data.frame(x = x[idx], y = y[idx], theta = theta[idx])
  }

  # final placer function
  function(lines, labels_data) {
    label_placer_simple(lines, labels_data, placer_fun)
  }
}


#' @rdname label_placer
#' @export
label_placer_middle <- function(rot_adjuster = angle_halfcircle_bottom()) {
  placer_fun <- function(line_data, ...) {
    out <- data.frame(x = numeric(0), y = numeric(0), theta = numeric(0))

    # It identifies each isoline subdivision. For an individual isoline the id column identifies the number of subdivisions.
    line_sections <- unique(line_data$id)

    # Then the label is printed at the middle of each isoline subdivision.
    for (i in 1:length(line_sections)) {
      x <- line_data$x[line_data$id == i]
      y <- line_data$y[line_data$id == i]

      middle_index <- as.integer(length(x) / 2)

      x_mid <- x[middle_index]
      y_mid <- y[middle_index]

      xtheta <- c(x[middle_index - 1], x[middle_index], x[middle_index + 1])
      ytheta <- c(y[middle_index - 1], y[middle_index], y[middle_index + 1])

      m <- cbind(xtheta - mean(xtheta), ytheta - mean(ytheta))
      v <- svd(m)$v

      out[i, ] <- list(x = x_mid, y = y_mid, theta = atan2(v[2], v[1]))
    }
    # standardize rotation angles for text labels
    out$theta <- rot_adjuster(out$theta)
    out
  }

  # final placer function
  function(lines, labels_data) {
    label_placer_simple(lines, labels_data, placer_fun)
  }
}

#' Standardize label angles
#'
#' Function factories that return functions to standardize rotation angles to specific angle ranges.
#'
#' `angle_halfcircle_bottom()` standardizes angles to (-pi/2, pi/2].
#'
#' `angle_halfcircle_right()` standardizes angles to (0, pi].
#'
#' `angle_fixed()` sets all angles to a fixed value (0 by default).
#'
#' `angle_identity()` does not modify any angles.
#' @param theta Fixed angle, in radians.
#' @export
angle_halfcircle_bottom <- function() {
  function(theta) {
    ifelse(
      theta <= -pi/2,
      theta + pi,
      ifelse(
        theta > pi/2,
        theta - pi,
        theta
      )
    )
  }
}

#' @rdname angle_halfcircle_bottom
#' @export
angle_halfcircle_right <- function() {
  function(theta) {
    ifelse(
      theta <= 0,
      theta + pi,
      ifelse(
        theta > pi,
        theta - pi,
        theta
      )
    )
  }
}

#' @rdname angle_halfcircle_bottom
#' @export
angle_fixed <- function(theta = 0) {
  force(theta)

  function(x) {
    rep_len(theta, length(x))
  }
}

#' @rdname angle_halfcircle_bottom
#' @export
angle_identity <- function() {
  function(x) {
    x
  }
}
