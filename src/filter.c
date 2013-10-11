/*******************************************************************************
 *  filter.c
 *  
 *  Created by Daniel Gatti on 3/14/11.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* a: 3D matrix of transition probilities for current DO generation.
 *             num states * num states * num SNPs.
 *  double* b: 3D matrix of emission log-probabilities. num states * num samples
 *              * num SNPs.
 *  double* prpred: 3D matrix of predictive log-probabilities. num states *
 *                  num samples * num SNPs.  We expect the initial values to
 *                  be in the first SNP slice.
 *  double* prfilt: 3D matrix of filtered log-probabilities. num states *
 *                  num samples * num SNPs.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through snps, then samples, then states should minimize cache-misses.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h> 
#include "addlog.h"


void filter_from_r(int* dims, double* a, double* b, double* prpred, 
		           double* prfilt) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st1 = 0; /* index for states */
  int st2 = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of samples in prpred & prfilt. */
  int num_samples = dims[1];
  int num_snps    = dims[2];
  register int outerindex = 0;  /* Index into 3D arrays. */
  register int innerindex = 0;  /* Index into 3D arrays. */
  int tr_mat_SNP = 0;           /* SNP index into trans. mat. */
  int pr_mat_slice = num_states * num_samples;  /* One SNP slice in prob. matrix. */
  int tr_mat_slice = num_states * num_states;   /* One SNP slice in trans. matrix. */
  double sums = -DBL_MAX;  /* Sums for each sample. */

  /* Initialize filtered probabilities by copying the predictive values at the
   * first SNP into the first SNP of the filtered array. */
  snp = 0;
  for(sam = 0; sam < num_samples; sam++) {
	/* Note that we have to initialize with this to make the addlog values correct.
	 * exp(0) is 1. exp(-Inf) = 0. */
	sums = -DBL_MAX;
	/* outerindex moves us to the current sample. */
	outerindex = sam * num_states;
    for(st1 = 0; st1 < num_states; st1++) {
      /* innerindex moves through states in one sample. */
	  innerindex = st1 + outerindex;
	  /* Numerator of Eqn. 11 in Churchill, 1989.  */
      prfilt[innerindex] = b[innerindex] + prpred[innerindex];
      /* Create the denominator for Eqn. 11 in Churchill, 1989. */
	  sums = addlog(sums, prfilt[innerindex]);
	} /* for(st1) */
	
	/* Normalize this sample by dividing by all states. */
	for(st1 = 0; st1 < num_states; st1++) {
	  /* Divide each state by the denominator for Eqn. 11 in Churchill, 1989. */
	  prfilt[st1 + outerindex] -= sums;
	} /* for(st1) */
  } /* for(sam) */
  
  /* Loop through each SNP, filling in the predictive array and then the 
   * filtered array. */
  for(snp = 1; snp < num_snps; snp++) {
	if(snp % 100 == 0)	R_CheckUserInterrupt();
    /* The transition matrix array has one less SNP than the probability
     * arrays.  So the index for the transition matrix is the SNP - 1. */
	/* tr_mat_SNP moves us to the current SNP in the transition matrix array. */
    tr_mat_SNP = (snp - 1) * tr_mat_slice;
    for(sam = 0; sam < num_samples; sam++) {
      sums = -DBL_MAX;
      /* outerindex moves us to the current sample and SNP in the probability 
       * arrays. */
	  outerindex = sam * num_states + snp * pr_mat_slice;
	  for(st1 = 0; st1 < num_states; st1++) {
		/* innerindex moves us to the current state. */ 
        innerindex = st1 + outerindex;
		prpred[innerindex] = -DBL_MAX;
		prfilt[innerindex] = -DBL_MAX;
		for(st2 = 0; st2 < num_states; st2++) {
          /* Eqn 10 in Churchill, 1989. */
		  prpred[innerindex] = addlog(prpred[innerindex], 
		                       a[st2 + st1 * num_states + tr_mat_SNP] +
							   prfilt[st2 + sam * num_states +
							   (snp - 1) * pr_mat_slice]);
        } /* for(st2) */
		/* Numerator in Eqn 11 in Churchill, 1989. */
        prfilt[innerindex] = b[innerindex] + prpred[innerindex];
        sums = addlog(sums, prfilt[innerindex]);
	  } /* for(st1) */

      /* Normalize this sample by dividing by all states. */
	  for(st1 = 0; st1 < num_states; st1++) {
        /* Divide by the denominator in Eqn 11 in Churchill, 1989. */
	    prfilt[st1 + outerindex] -= sums;
	  } /* for(st1) */
	} /* for(sam) */
  } /* for(snp) */
} /* filter_from_r() */

