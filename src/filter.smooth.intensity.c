/*******************************************************************************
 *  filter.smooth.intensity.c
 *  
 *  Created by Daniel Gatti on 3/14/11.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* a: 3D matrix of transition probilities for current DO generation.
 *             num states * num states * num SNPs.
 *  double* b: 3D matrix of emission log-probabilities. num states * num samples
 *              * num SNPs.
 *  double* prsmth: 3D matrix of smoothed log-probabilities. num states *
 *                  num samples * num SNPs.  We expect the initial values to
 *                  be in the first SNP slice.
 *  double* init: Vector of initial probabilities of length num_states.
 *  double* loglik: double in which to return the log-likelihood.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through snps, then samples, then states should minimize cache-misses.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h> 
#include "addlog.h"

void filter_smooth_intensity(int* dims, double* a, double* b, double* prsmth,
                             double* init, double* loglik) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st1 = 0; /* index for states */
  int st2 = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of samples in prpred & prfilt. */
  int num_samples = dims[1];
  int num_snps    = dims[2];
  int pr_snp_index = 0;  /* SNP index into 3D arrays. */
  int pr_sam_index = 0;  /* Sample index into 3D arrays. */
  int pr_st_index = 0;   /* State index into 3D arrays. */
  int tr_snp_index = 0;    /* SNP index into trans. mat. */
  int pr_mat_slice = num_states * num_samples;  /* One SNP slice in prob. matrix. */
  int tr_mat_slice = num_states * num_states;   /* One SNP slice in trans. matrix. */
  double sums = -DBL_MAX;  /* Sums for each sample. */
  double* prpred = (double*)malloc(num_states * num_samples * num_snps * sizeof(double));
  double* prfilt = (double*)malloc(num_states * num_samples * num_snps * sizeof(double));

  /*** Filtering (or forward) ***/
  
  /* Initialize filtered probabilities by copying the predictive values at the
   * first SNP into the first SNP of the filtered array. */
  loglik[0] = -DBL_MAX;
  snp = 0;

  for(sam = 0; sam < num_samples; sam++) {
    /* Note that we have to initialize with this to make the addlog values correct.
     * exp(0) is 1. exp(-Inf) = 0. */
    sums = -DBL_MAX;
    /* pr_sam_index moves us to the current sample. */
    pr_sam_index = sam * num_states;
    for(st1 = 0; st1 < num_states; st1++) {
      /* pr_st_index moves through states in one sample. */
      pr_st_index = st1 + pr_sam_index;
      prpred[pr_st_index] = init[st1];
      /* Numerator of Eqn. 11 in Churchill, 1989.  */
      prfilt[pr_st_index] = b[pr_st_index] + prpred[pr_st_index];
      loglik[0] = addlog(loglik[0], prfilt[pr_st_index]);
      /* Create the denominator for Eqn. 11 in Churchill, 1989. */
      sums = addlog(sums, prfilt[pr_st_index]);
    } /* for(st1) */
	
    /* Normalize this sample by dividing by all states. */
    for(st1 = 0; st1 < num_states; st1++) {
      /* Divide each state by the denominator for Eqn. 11 in Churchill, 1989. */
      prfilt[st1 + pr_sam_index] -= sums;
    } /* for(st1) */
  } /* for(sam) */
     
  /* Loop through each SNP, filling in the predictive array and then the 
   * filtered array. */
  for(snp = 1; snp < num_snps; snp++) {
    if(snp % 100 == 0)	R_CheckUserInterrupt();
	
    /* The transition matrix array has one less SNP than the probability
     * arrays.  So the index for the transition matrix is the SNP - 1. */
	/* tr_snp_index moves us to the current SNP in the transition matrix array. */
    tr_snp_index = (snp - 1) * tr_mat_slice;
    pr_snp_index = snp * pr_mat_slice;
    for(sam = 0; sam < num_samples; sam++) {
      sums = -DBL_MAX;
      /* pr_sam_index moves us to the current sample and SNP in the probability 
       * arrays. */
      pr_sam_index = sam * num_states + pr_snp_index;
      for(st1 = 0; st1 < num_states; st1++) {
        /* pr_st_index moves us to the current state. */ 
        pr_st_index = st1 + pr_sam_index;
        prpred[pr_st_index] = -DBL_MAX;
        prfilt[pr_st_index] = -DBL_MAX;
        for(st2 = 0; st2 < num_states; st2++) {
          /* Eqn 10 in Churchill, 1989. */
          prpred[pr_st_index] = addlog(prpred[pr_st_index], 
                                a[st2 + st1 * num_states + tr_snp_index] +
                                prfilt[st2 + sam * num_states +
                                (snp - 1) * pr_mat_slice]);
        } /* for(st2) */
        /* Numerator in Eqn 11 in Churchill, 1989. */
        prfilt[pr_st_index] = b[pr_st_index] + prpred[pr_st_index];
        loglik[0] = addlog(loglik[0], prfilt[pr_st_index]);
        sums = addlog(sums, prfilt[pr_st_index]);
      } /* for(st1) */

      /* Normalize this sample by dividing by all states. */
      for(st1 = 0; st1 < num_states; st1++) {
        /* Divide by the denominator in Eqn 11 in Churchill, 1989. */
        prfilt[st1 + pr_sam_index] -= sums;
      } /* for(st1) */
    } /* for(sam) */
  } /* for(snp) */

  /*** Smoothing (or backward) ***/
  
  /* Initialize filtered probabilities by copying the filtered values at the
   * last SNP into the last SNP of the smoothed array. */
  snp = num_snps - 1;
  for(sam = 0; sam < num_samples; sam++) {
    /* pr_sam_index moves us to the current SNP and sample in the probability
     * arrays. */
    pr_sam_index = sam * num_states + snp * pr_mat_slice;
    for(st1 = 0; st1 < num_states; st1++) {
      /* pr_st_index moves us to the correct state. */
      pr_st_index = st1 + pr_sam_index;
      prsmth[pr_st_index] = prfilt[pr_st_index];
    } /* for(st1) */
  } /* for(sam) */

  /* Loop backward through all of the SNPs and update the smoothed probabilities. */
  for(snp = num_snps - 2; snp >= 0; snp--) {
    if(snp % 100 == 0) R_CheckUserInterrupt();
    /* tr_snp_index moves us to the correct SNP in the transition matrix array. */
    tr_snp_index = snp * tr_mat_slice;
    pr_snp_index = (snp + 1) * pr_mat_slice;
    for(sam = 0; sam < num_samples; sam++) {
      /* pr_sam_index moves us to the successive SNP and sample in the probability
       * arrays. */
      pr_sam_index = sam * num_states + pr_snp_index;
      for(st1 = 0; st1 < num_states; st1++) {
        sums = -DBL_MAX;
        for(st2 = 0; st2 < num_states; st2++) {
          /* pr_st_index moves us to the current state. */
	  pr_st_index = st2 + pr_sam_index;
	  /* Left summation of Eqn. 13 in Churchill, 1989 */
	  sums = addlog(sums, a[st1 + st2 * num_states + tr_snp_index] +
	         prsmth[pr_st_index] - prpred[pr_st_index]);
	} /* for(st2)*/
	/* pr_st_index moves us to the current SNP and sample in the probability
	 * arrays. */
	pr_st_index = st1 + sam * num_states + snp * pr_mat_slice;
	/* Multiply by prfilt as in Eqn. 13 of Churchill, 1989. */
	prsmth[pr_st_index] = prfilt[pr_st_index] + sums;
      } /* for(st1) */
    } /* for(sam) */
  } /* for(snp) */
  
  free(prpred);
  free(prfilt);

} /* filter_smooth_intensity() */


