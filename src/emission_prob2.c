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

void emission_prob2(int* dims, double* x, double* y, double* xmeans,
                   double* ymeans, double* xvars, double* yvars, double* covars, 
                   double* probs) {

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
  int x_index = 0;           /* Index into x & y arrays. */
  int x_snp_index = 0;  /* Index into x & y arrays. */
  double x_diff = 0.0;  /* Difference between x & x mean. */
  double y_diff = 0.0;  /* Difference between y & y mean. */
  double covar_x_diff_y_diff = 0.0; /* Product of covar * x_diff * y_diff. */
  double det    = 0.0;  /* Determinant of cluster covariance matrix.*/
  double sample_sum = 0.0; /* Sum of probs for each sample. */
  double default_prob = log(1.0 / (double)num_states); /* Uniform prob for
                                                          missing data.  */
  double log2pi = log(2.0 * M_PI);  /* log(2 * PI) */

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

    /* x_snp_index moves us to the correct SNP slice in x & y. */
    x_snp_index = snp * num_samples;
 
    /* Calculate emission probabilities for all samples. */
    for(sam = 0; sam < num_samples; sam++) {
      
      /* Zero out the sample sum accumulator. */
      sample_sum = -DBL_MAX;
      
      /* prob_sample_index moves us to the sample column in probs. */
      prob_sample_index = prob_snp_index + sam * num_states;

      /* x_index moves us to the current sample in the x & y matrices. */
      x_index = sam + x_snp_index;

	 /* Set the intensity < 0.0 to indicate missing data. */
	 if((x[x_index] >= 0.0) && (y[x_index] >= 0.0)) {

        /* Calculate the emission probabilities for each state. */
        for(st = 0; st < num_states; st++) {

          /* state_index moves us to the current state in the mean & var arrays. */
          state_index = st + snp_index;
          prob_state_index = st + prob_sample_index;
        
          /* Calculate x and y difference from the state mean. */
          x_diff = x[x_index] - xmeans[state_index];
          y_diff = y[x_index] - ymeans[state_index];
          covar_x_diff_y_diff = covars[state_index] * x_diff * y_diff;

          /* Calculate the emission probability from the bivariate Gaussian
           * distribution. */
          /* Determinant of covariance matrix. */
          det = xvars[state_index] * yvars[state_index] - covars[state_index] * covars[state_index];

          /* The inverse of a 2x2 matrix |ab| is 1/det | d -b|
           *                             |cd|          |-c  a| */
          /* This is a highly factored version of the bivariate Gaussian density. */
          probs[prob_state_index] = log(sqrt(det)) - log2pi - (0.5 / det * 
                (yvars[state_index] * x_diff * x_diff + xvars[state_index] * y_diff * y_diff - 
                2 * covar_x_diff_y_diff));

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
} /* emission_prob2() */
