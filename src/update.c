/*******************************************************************************
 *  update.c
 *  
 *  Created by Daniel Gatti on Aug. 16, 2011.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* t: matrix of theta values. num samples * num SNPs.
 *  double* r: matrix of rho values. num samples * num SNPs.
 *  double* tmeans: 2D matrix of theta and rho means. num states * num SNPs.
 *  double* rmeans: 2D matrix of theta and rho means. num states * num SNPs.
 *  double* tvars: 2D matrix of theta variances. num states * num SNPs.
 *  double* rvars: 2D matrix of rho variances. num states * num SNPs.
 *  double* prsmth: 3D matrix of smoothed log-probabilities. num states *
 *                  num samples * num SNPs.
 *  double* foundertmeans: 2D matrix of founder theta mean intensities.
 *                         num_states * num_snps.
 *  double* founderrmeans: 2D matrix of founder rho mean intensities.
 *                         num_states * num_snps.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through states, then samples, then SNPs should minimize cache-misses.
 * At each SNP, the variance is pooled among all states and samples.
 * We perform the calculations on a non-log scale.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include <math.h>

void update_intensity(int* dims, double* t, double* r, double* tmeans, 
                   double* rmeans, double* tvars, double* rvars,
                   double* prsmth, double* foundertmeans, double* founderrmeans) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st  = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of states in prsmth. */
  int num_samples = dims[1]; /* Total number of samples in prsmth. */
  int num_snps    = dims[2]; /* Total number of SNPs in prsmth. */
  int snp_index = 0;    /* SNP index for mean and covar arrays. */
  int state_index = 0;  /* State index for mean and covar arrays. */
  int t_snp_index = 0;  /* Index into theta & rho arrays. */
  int t_index = 0;      /* Index into theta & rho arrays. */
  double t_diff = 0.0;  /* Difference between theta & theta mean. */
  double r_diff = 0.0;  /* Difference between rho & rho mean. */
  int prsmth_slice = num_states * num_samples; /* One SNP slice in prsmth. */
  int prsmth_snp_index = 0; /* SNP index into prsmth. */
  double prsmth_sums[num_states];  /* prsmth sums for each state. */
  double exp_prsmth[num_samples]; /* exp(prsmth) for the current sample. */
  int pad = 2.0;  /* The number of founder samples that we use to pad
                     the theta and rho locations. */ 
  
  /* Update the state means and variances. 
   * We use all samples to update the means and variances.
   * This is because if we have founders and F1s, we want their values to 
   * anchor the means. */
  for(snp = 0; snp < num_snps; snp++) {
    if(snp % 100 == 0) R_CheckUserInterrupt();

    /* Update the state means. 
	 * snp_index moves us to the current SNP in the mean arrays. */
	snp_index   = snp * num_states;
	t_snp_index = snp * num_samples;
	prsmth_snp_index = snp * prsmth_slice;

	/* prsmth_sums will hold the sum of prsmth values for each state across all 
	 * samples at the current SNP. */
    for(st = 0; st < num_states; st++) {
      /* Pad the counts w/ the founder means. */
      prsmth_sums[st] = pad;
      
	  /* state_index moves us to the current state in the mean arrays. */
      state_index = st + snp_index;
      tmeans[state_index] = pad * foundertmeans[state_index];
      rmeans[state_index] = pad * founderrmeans[state_index];

	  /* Use all samples to calculate the new state means. */
      for(sam = 0; sam < num_samples; sam++) {
        /* t_index moves us to the current sample in the theta & rho matrices. */
	    t_index = sam + t_snp_index;
		/* Only use this data point if it is > 0.0. Negative values are missing data. */
		if((t[t_index] >= 0.0) & (r[t_index] >= 0.0)) {
          exp_prsmth[sam] = exp(prsmth[st + sam * num_states + prsmth_snp_index]);
	      /* Update the theta & rho means. */
          tmeans[state_index] += t[t_index] * exp_prsmth[sam];
          rmeans[state_index] += r[t_index] * exp_prsmth[sam];
          /* Add the current prsmth values to the sum for all samples at this
           * state. */
		  prsmth_sums[st] += exp_prsmth[sam];
        } /* if(t[t_index] > 0.0 & r[t_index] > 0.0) */
      } /* for(sam) */

      /* Divide by the sum of the prsmth at this state. */
      tmeans[state_index] /= prsmth_sums[st];
      rmeans[state_index] /= prsmth_sums[st];

      /* Using the new means, update the state covariances. */
      tvars[state_index]  = 0.0; 
      rvars[state_index]  = 0.0; 
      prsmth_sums[st] -= pad;
      
      /* If there are not enough samples in the current state,
       * make sure the the denominator is not too small. */
      if(prsmth_sums[st] < 1e-16) {
    	  prsmth_sums[st] = 1e-16;
      } /* if(prsmth_sums[st] < 1e-16) */

      /* Use all samples to calculate the new state covariances. */
      for(sam = 0; sam < num_samples; sam++) {
        /* t_index moves us to the current sample in the theta & rho matrices. */
	    t_index = sam + t_snp_index;
		/* Only use this data point if it is > 0.0. Negative values are missing data. */
		if((t[t_index] >= 0.0) & (r[t_index] >= 0.0)) {	  	
  		  t_diff = t[t_index] - tmeans[state_index];
  		  tvars[state_index]  += t_diff * t_diff * exp_prsmth[sam];

	  	  r_diff = r[t_index] - rmeans[state_index];
	  	  rvars[state_index]  += r_diff * r_diff * exp_prsmth[sam];
        } /* if(t[t_index] > 0.0 & r[t_index] > 0.0) */

        /* Zero out the exp(prsmth) for this state. */
        exp_prsmth[sam] = 0.0;
      } /* for(sam) */

      /* Divide by sum of all probs at this state.*/
      state_index = st + snp_index;
      tvars[state_index]  /= prsmth_sums[st];
      rvars[state_index]  /= prsmth_sums[st];

    } /* for(st) */
  } /* for(snp) */
} /* update_intensity() */
