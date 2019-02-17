# convenience operator, don't want to import rlang just for this
"%||%" <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}
