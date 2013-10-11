/**********************************************************************
 * DOstep.h
 *
 * function declarations ... for calculating transition probabilities
 * for the diversity outcross
 *
 * Karl W Broman
 * first written 23 Aug 2011
 * last modified 23 Aug 2011
 *
 **********************************************************************/

double ri4way_auto_hapAA(double r, int k);
double ri4way_femX_hapAA(double r, int k);
double ri4way_femX_hapCC(double r, int k);
double ri4way_malX_hapAA(double r, int k);
double ri4way_malX_hapCC(double r, int k);

double DOrec_auto(double r, int s, int n_k, int *k, double *alpha_k);
double DOrec_femX_s1(double r, int n_k, int *k, double *alpha_k);
double DOrec_malX_s1(double r, int n_k, int *k, double *alpha_k);
double DOrec_femX(double r, int s, int n_k, int *k, double *alpha_k);
double DOrec_malX(double r, int s, int n_k, int *k, double *alpha_k);

double DOstep_auto(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k);
double DOstep_femX(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k);
double DOstep_malX(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k);

/* end of DOstep.h */
