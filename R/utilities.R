# convenience operator, don't want to import rlang just for this
"%||%" <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

# evaluates all arguments
# (simpler than forcing each argument individually)
force_all <- function(...) list(...)

rethrow_interrupt <- function() {
  interrupt <- structure(list(), class = c("interrupt", "condition"))
  signalCondition(interrupt)
  invokeRestart("abort")
}
