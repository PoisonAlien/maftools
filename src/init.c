#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// External function declarations
extern SEXP bamrc_c(SEXP bam_path, SEXP chr_vec, SEXP pos_vec, SEXP mapq, SEXP flag, SEXP fasta_path, SEXP include_idxstats, SEXP verbose, SEXP hotspot_data);

// Registration table
static const R_CallMethodDef callMethods[] = {
    {"bamrc_c", (DL_FUNC) &bamrc_c, 9},
    {NULL, NULL, 0}
};

// Initialize the dynamic library
void R_init_maftools(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}