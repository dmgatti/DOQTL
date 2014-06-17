/**********************************************************************
 * DOstep_Rwrappers.h
 * 
 * wrappers for calling various functions related to calculate 
 * transitions probabilities for DO from R
 *
 **********************************************************************/

void R_ri4way_auto_hapAA(double *r, int *k, double *result);
void R_ri4way_femX_hapAA(double *r, int *k, double *result);
void R_ri4way_femX_hapCC(double *r, int *k, double *result);
void R_ri4way_malX_hapAA(double *r, int *k, double *result);
void R_ri4way_malX_hapCC(double *r, int *k, double *result);

void R_DOrec_auto(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result);
void R_DOrec_femX(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result);
void R_DOrec_malX(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result);

void R_DOstep_auto(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result);
void R_DOstep_femX(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result);
void R_DOstep_malX(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result);

/* end of DOstep_Rwrappers.h */
