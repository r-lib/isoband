#' Set up a label placement strategy
#'
#' The position and rotation of labels is calculated based on the line
#' data for individual isolines.
#' @param placement String consisting of any combination of the letters
#' "t", "r", "b", "l" indicating the placement of labels at the top,
#' to the right, at the bottom, to the left of the isoline.
#' @param rot_adjuster Function that standardizes the rotation angles of the labels.
#' @param n Size of the point neighborhood over which the rotation angle should be
#' calculated.
#' @export
label_placer_minmax <- function(placement = "tb", rot_adjuster = angle_halfcircle_bottom, n = 2) {
  force_all(placement, rot_adjuster, n)

  function(line_data) {
    # find location for labels
    idx <- na.omit(
      c(
        which( # placement "top"
          isTRUE(grepl("t", placement, fixed = TRUE)) & line_data$y == min(line_data$y)
        )[1],
        which( # placement "bottom"
          isTRUE(grepl("b", placement, fixed = TRUE)) & line_data$y == max(line_data$y)
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


#' Standardize label angles
#'
#' Functions to standardize rotation angles to specific angle ranges.
#'
#' `angle_halfcircle_bottom()` standardizes angles to (-pi/2, pi/2].
#' `angle_halfcircle_right()` standardizes angles to (0, pi].
#' `angle_halfcircle_zero()` sets all angles to 0.
#' `angle_identity()` does not modify any angles.
#' @param theta Numeric vector of angles in radians.
#' @export
angle_halfcircle_bottom <- function(theta) {
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

#' @rdname angle_halfcircle_bottom
#' @export
angle_halfcircle_right <- function(theta) {
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

#' @rdname angle_halfcircle_bottom
#' @export
angle_zero <- function(theta) {
  rep_len(0, length(theta))
}

#' @rdname angle_halfcircle_bottom
#' @export
angle_identity <- function(theta) {
  theta
}
