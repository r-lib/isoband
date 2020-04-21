#ifndef UTILS_H
#define UTILS_H

#define BEGIN_CPP try {

#define END_CPP                                                                \
}                                                                              \
catch (std::exception & e) {                                                   \
  Rf_error("C++ exception: %s", e.what());                                     \
}

#endif
