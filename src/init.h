void emission_prob(int* dims, double* theta, double* rho, double* thetameans, 
                          double* rhomeans, double* thetavars, double* rhovars, double* probs);

void filter_smooth_intensity(int* dims, double* a, double* b, double* prsmth,
                             double* init, double* loglik);
							 
void filter_smooth_allele(int* dims, int* geno, double* a, double* b, double* prsmth,
                          double* init, double* loglik);

void kinship(int* dims, double* probs, double* K);

void update_intensity(int* dims, double* t, double* r, double* tmeans, 
                   double* rmeans, double* tvars, double* rvars,
                   double* prsmth, double* foundertmeans, double* founderrmeans);

void update_alleles(int* dims,int* geno, double* b, double* pseudo, 
                           double* prsmth);

void viterbi(int* dims, double* a, double* b, int* v_path, double* v_prob);

void DO_autosome_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k, 
		                     int* num_recomb, double* recomb);

void DO_femaleX_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k,
                            int* num_recomb, double* recomb);

void DO_maleX_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k,
	                      int* num_recomb, double* recomb);
