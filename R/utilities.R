# evaluates all arguments
# (simpler than forcing each argument individually)
force_all <- function(...) list(...)

rethrow_interrupt <- function() {
  interrupt <- structure(list(), class = c("interrupt", "condition"))
  signalCondition(interrupt)
  invokeRestart("abort")
}
