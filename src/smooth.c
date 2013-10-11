/*******************************************************************************
 *  smooth.c
 *  
 *  Created by Daniel Gatti on 3/14/11.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of SNPs, number of states,
 *             number of samples.
 *  double* a: matrix of transition probilities. num states * num states.
 *  double* prpred: 3D matrix of predictive log-probabilities. num states *
 *                  num samples * num SNPs.
 *  double* prfilt: 3D matrix of filtered log-probabilities. num states *
 *                  num samples * num SNPs.
 *  double* prsmth: 3D matrix of smoothed log-probabilities. num states *
 *                  num samples * num SNPs.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through snps, then samples, then states should minimize cache-misses.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "addlog.h"

void smooth_from_r(int* dims, double* a, double* prpred, double* prfilt,
                   double* prsmth) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st1 = 0; /* index for states */
  int st2 = 0; /* index for states */
  int num_states  = dims[0]; /* Number of states in data. */
  int num_samples = dims[1]; /* Total number of samples in data. */
  int num_snps    = dims[2]; /* Number of SNPs in data. */
  register int outerindex = 0;  /* Index into 3D arrays. */
  register int innerindex = 0;  /* Index into 3D arrays. */
  int tr_mat_SNP = 0;           /* SNP index into trans. mat. */
  int pr_mat_slice = num_samples * num_states; /* The size of one slice in the 3D array. */
  int tr_mat_slice = num_states *  num_states; /* One SNP slice in trans. matrix. */
  double sums = -DBL_MAX;  /* Sums for each state. */
  
  /* Initialize filtered probabilities by copying the filtered values at the
   * last SNP into the last SNP of the smoothed array. */
  snp = num_snps - 1;
  for(sam = 0; sam < num_samples; sam++) {
    /* outerindex moves us to the current SNP and sample in the probability
     * arrays. */
    outerindex = sam * num_states + snp * pr_mat_slice;
    for(st1 = 0; st1 < num_states; st1++) {
      /* innerindex moves us to the correct state. */
	  innerindex = st1 + outerindex;
      prsmth[innerindex] = prfilt[innerindex];
	} /* for(st1) */
  } /* for(sam) */

  /* Loop backward through all of the SNPs and update the smoothed probabilities. */
  for(snp = num_snps - 2; snp >= 0; snp--) {
	if(snp % 100 == 0) R_CheckUserInterrupt();
    /* tr_mat_SNP moves us to the correct SNP in the transition matrix array. */
	tr_mat_SNP = snp * tr_mat_slice;
    for(sam = 0; sam < num_samples; sam++) {
      /* outerindex moves us to the succesive SNP and sample in the probability
       * arrays. */
	  outerindex = sam * num_states + (snp + 1) * pr_mat_slice;
	  for(st1 = 0; st1 < num_states; st1++) {
        sums = -DBL_MAX;
		for(st2 = 0; st2 < num_states; st2++) {
          /* innerindex moves us to the current state. */
		  innerindex = st2 + outerindex;
		  /* Left summation of Eqn. 13 in Churchill, 1989 */
		  sums = addlog(sums, a[st1 + st2 * num_states + tr_mat_SNP] +
				        prsmth[innerindex] - prpred[innerindex]);
		} /* for(st2)*/
		/* innerindex moves us to the current SNP and sample in the probability
		 * arrays. */
		innerindex = st1 + sam * num_states + snp * pr_mat_slice;
		/* Multiply by prfilt as in Eqn. 13 of Churchill, 1989. */
		prsmth[innerindex] = prfilt[innerindex] + sums;
	  } /* for(st1) */
	} /* for(sam) */
  } /* for(snp) */
} /* smooth_from_r() */

