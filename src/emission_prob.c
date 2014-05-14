/*******************************************************************************
 *  emission.prob.c
 *
 *  Created by Daniel Gatti on Dec. 5, 2012.
 *  Copyright 2012 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* theta: matrix of theta values. num states * num SNPs.
 *  double* rho: matrix of rho values. num states * num SNPs.
 *  double* thetameans: matrix of theta and rho means. num states * num SNPs.
 *  double* rhomeans: matrix of theta and rho means. num states * num SNPs.
 *  double* thetavars:  matrix of theta variances. num states * num SNPs * 2.
 *  double* rhovars:  matrix of rho variances. num states * num SNPs * 2.
 *  double* probs:  3D matrix of emission probabilities that are returned.
 *                  num.states * num.samples * num.SNPs.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through states, then samples, then SNPs should minimize cache-misses.
 * At each SNP, the variance is pooled among all states and samples.
 * We perform the calculations on a non-log scale.
 */
#include <R.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h> 
#include "addlog.h"

void emission_prob(int* dims, double* theta, double* rho, double* thetameans,
                   double* rhomeans, double* thetavars, double* rhovars, double* probs) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st  = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of states in prsmth. */
  int num_samples = dims[1]; /* Total number of samples in prsmth. */
  int num_snps    = dims[2]; /* Total number of SNPs in prsmth. */
  int snp_index = 0;    /* SNP index for mean and covar arrays. */
  int state_index = 0;  /* State index for mean and covar arrays. */
  int prob_slice = num_states * num_samples;   /* The size of one SNP slice in 
                                               probs.*/
  int prob_snp_index = 0;    /* Index into probs array. */
  int prob_sample_index = 0; /* Index into probs array. */
  int prob_state_index = 0;  /* Index into probs array. */
  int theta_index = 0;         /* Index into theta & y arrays. */
  int theta_snp_index = 0;  /* Index into theta & y arrays. */
  double theta_diff = 0.0;  /* Difference between theta & theta mean. */
  double rho_diff = 0.0;  /* Difference between rho & rho mean. */
  double sample_sum = 0.0; /* Sum of probs for each sample. */
  double default_prob = log(1.0 / (double)num_states); /* Uniform prob for
                                                          missing data.  */

  /* Update the state means and variances.
   * We use all samples to update the means and variances.
   * This is because if we have founders and F1s, we want their values to
   * anchor the means. */
  for(snp = 0; snp < num_snps; snp++) {  
    /* Update the state means. */
	/* snp_index moves us to the current SNP in the mean arrays. */
	snp_index = snp * num_states;
	/* prob_snp_index moves us to the correct SNP slice in probs. */
	prob_snp_index = snp * prob_slice;
	/* theta_snp_index moves us to the correct SNP slice in theta & rho. */
    theta_snp_index = snp * num_samples;
 
	/* Calculate emission probabilities for all samples. */
    for(sam = 0; sam < num_samples; sam++) {
      
      /* Zero out the sample sum accumulator. */
      sample_sum = -DBL_MAX;
      
      /* prob_sample_index moves us to the sample column in probs. */
      prob_sample_index = prob_snp_index + sam * num_states;

      /* theta_index moves us to the current sample in the theta & rho matrices. */
      theta_index = sam + theta_snp_index;

	  /* Set the intensity < 0.0 to indicate missing data. */
	  if((theta[theta_index] >= 0.0) && (rho[theta_index] >= 0.0)) {

        /* Calculate the emission probabilities for each state. */
        for(st = 0; st < num_states; st++) {

          /* state_index moves us to the current state in the mean & var arrays. */
          state_index = st + snp_index;
          prob_state_index = st + prob_sample_index;
        
	      /* Calculate x and rho difference from the state mean. */
          theta_diff = theta[theta_index] - thetameans[state_index];
          rho_diff   = rho[theta_index]   - rhomeans[state_index];

          /* Calculate the emission probability. */
          probs[prob_state_index] = 
        		log(0.5 / M_PI / sqrt(thetavars[state_index] * rhovars[state_index])) -
        		0.5 * (((theta_diff * theta_diff) / thetavars[state_index]) + 
        		       ((rho_diff * rho_diff) / rhovars[state_index]));
          sample_sum = addlog(sample_sum, probs[prob_state_index]);
        } /* for(st) */
      } else {
	    /* We have no intensity value. Place a uniform probability in all states. */
        for(st = 0; st < num_states; st++) {

          prob_state_index = st + prob_sample_index;
 
          /* Calculate the emission probability. */
          probs[prob_state_index] = default_prob;
          sample_sum = addlog(sample_sum, probs[prob_state_index]);
        } /* for(st) */
	  } /* else */
		
      /* Divide by the sum of all states so that the probabilities sum to 1. */
      for(st = 0; st < num_states; st++) {
        probs[st + prob_sample_index] -= sample_sum;
      } /* for(st) */
      
    } /* for(sam) */      

  } /* for(snp) */  
} /* emission_prob() */
