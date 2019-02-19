# Calculate the label position for one set of isolines (one level). Used by
# `isolines_grob()` to place labels.
#
# @param line_data line segments, specified as list of x, y, id
# @param break_id character vector specifying the break identifier
# @param label character vector specifying the label that will be printed
#   instead of the break identifier
# @param index index into original list of all isolines
# @param break_index index into original list of breaks
#
# The parameters `break_id`, `label`, `index`, and `break_index` are provided
# simply so they can be added to the resulting data frame holding
# label positions
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
    out[i, ] <- label_position_impl(data, idx[i], n)
  }

  out
}

# function that does all the work for label_positions().
# requires a single index idx
label_position_impl <- function(data, idx, n = 2) {
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
