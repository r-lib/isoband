#ifndef UTILS_H
#define UTILS_H

// Turn off non-exported functionality
#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

// Define a C++ try-catch macro to guard C++ calls
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

class CollectorList {
  SEXP data_;
  R_xlen_t n_;

public:
  CollectorList(R_xlen_t size = 1) : n_(0) {
    data_ = Rf_allocVector(VECSXP, size);
    R_PreserveObject(data_);
  }

  void push_back(SEXP x) {
    if (Rf_xlength(data_) == n_) {
      R_ReleaseObject(data_);
      data_ = Rf_lengthgets(data_, n_ * 2);
      R_PreserveObject(data_);
    }
    SET_VECTOR_ELT(data_, n_++, x);
  }

  operator SEXP() {
    if (Rf_xlength(data_) != n_) {
      SETLENGTH(data_, n_);
    }
    return data_;
  }

  ~CollectorList() { R_ReleaseObject(data_); }
};

#endif
