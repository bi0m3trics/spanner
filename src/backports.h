#ifndef SPANNER_BACKPORTS_H
#define SPANNER_BACKPORTS_H

#include <R.h>
#include <Rinternals.h>

// Backports for R < 4.5.0 compatibility
// Based on Writing R Extensions manual section 6.22

#if R_VERSION < R_Version(4, 5, 0)

// CLEAR_ATTRIB was added in R 4.5.0
#ifndef CLEAR_ATTRIB
static inline void CLEAR_ATTRIB(SEXP x) {
    SET_ATTRIB(x, R_NilValue);
    SET_OBJECT(x, 0);
    UNSET_S4_OBJECT(x);
}
#endif

// ANY_ATTRIB was added in R 4.5.0
#ifndef ANY_ATTRIB
static inline int ANY_ATTRIB(SEXP x) {
    return ATTRIB(x) != R_NilValue;
}
#endif

#endif // R_VERSION < R_Version(4, 5, 0)

#if R_VERSION < R_Version(4, 6, 0)

// R_class was added in R 4.6.0
#ifndef R_class
#define R_class(x) R_data_class(x, FALSE)
#endif

// DATAPTR_RW was added in R 4.6.0
#ifndef DATAPTR_RW
#define DATAPTR_RW(x) DATAPTR(x)
#endif

// Resizable vector functions were added in R 4.6.0
#ifndef R_allocResizableVector
static inline SEXP R_allocResizableVector(SEXPTYPE type, R_xlen_t maxlen) {
    SEXP ret = Rf_allocVector(type, maxlen);
    return ret;
}
#endif

#ifndef R_isResizable
static inline bool R_isResizable(SEXP x) {
    return false;  // Vectors are not resizable in R < 4.6.0
}
#endif

#ifndef R_maxLength
static inline R_xlen_t R_maxLength(SEXP x) {
    return Rf_xlength(x);
}
#endif

#ifndef R_resizeVector
static inline void R_resizeVector(SEXP x, R_xlen_t newlen) {
    // In older R versions, we cannot safely resize vectors
    // This is a placeholder that should not be called
    Rf_error("Vector resizing not supported in R < 4.6.0");
}
#endif

#ifndef R_duplicateAsResizable
static inline SEXP R_duplicateAsResizable(SEXP x) {
    return Rf_duplicate(x);
}
#endif

#endif // R_VERSION < R_Version(4, 6, 0)

#endif // SPANNER_BACKPORTS_H
