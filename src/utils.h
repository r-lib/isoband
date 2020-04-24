#ifndef UTILS_H
#define UTILS_H

#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

#define BEGIN_CPP try {

#define END_CPP                                                                \
}                                                                              \
catch (std::exception & e) {                                                   \
  Rf_error("C++ exception: %s", e.what());                                     \
}

static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
static inline bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

[[ noreturn ]] static inline void longjump_interrupt() {
  SEXP env = PROTECT(Rf_findVarInFrame(R_NamespaceRegistry, Rf_install("isoband")));
  if (env == R_UnboundValue) {
    Rf_error("isoband namespace could not be found");
  }
  SEXP call = PROTECT(Rf_lang1(Rf_install("rethrow_interrupt")));
  Rf_eval(call, env);
  Rf_error("Interrupt failed to rethrow");
}

#endif
