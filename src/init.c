#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "init.h"

static const R_CMethodDef cMethods[]  = {
    {"emission_prob", (DL_FUNC) &emission_prob, 8},
    {"emission_prob2", (DL_FUNC) &emission_prob2, 9},
    {"filter_smooth_intensity", (DL_FUNC) &filter_smooth_intensity, 6},
    {"filter_smooth_allele", (DL_FUNC) &filter_smooth_allele, 7},
    {"kinship", (DL_FUNC) &kinship, 3},
    {"update_intensity", (DL_FUNC) &update_intensity, 10},
    {"update_intensity2", (DL_FUNC) &update_intensity2, 11},
    {"update_alleles", (DL_FUNC) &update_alleles, 5},
    {"viterbi", (DL_FUNC) &viterbi, 5},
    {"DO_autosome_recomb_freq", (DL_FUNC) &DO_autosome_recomb_freq, 7},
    {"DO_maleX_recomb_freq", (DL_FUNC) &DO_maleX_recomb_freq, 7},
    {"DO_femaleX_recomb_freq", (DL_FUNC) &DO_femaleX_recomb_freq, 7},
    {NULL, NULL, 0}
};


void attribute_visible R_init_DOQTL(DllInfo *info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

