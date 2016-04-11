/*******************************************************************************
 *  kinship.c
 *  
 *  Created by Daniel Gatti on 8/31/2013.
 *  Copyright 2013 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* probs: 3D matrix of founder probabilities. num samples * 
 *                 num founders * num SNPs.
 *  double* K: 2D matrix of kinships. num samples * num samples.
 * We use the relation cos(theta) = A.B / (||A||*||B||)
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through snps, then samples, then states should minimize cache-misses.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h> 
#include <string.h>
#include "addlog.h" 

void kinship(int* dims, double* probs, double* K) {

  int snp = 0; /* index for SNPs */
  int sam1 = 0; /* index for samples */
  int sam2 = 0; /* index for samples */
  int f1 = 0; /* index for founders */
  int num_samples  = dims[0]; /* Total number of samples in probs * K. */
  int num_founders = dims[1]; 
  int num_snps     = dims[2];
  int pr_mat_slice = num_founders * num_samples;  /* One SNP slice in prob. matrix. */
  int pr_mat_snp_index = 0; /* SNP index into probs matrix. */
  int pr_mat_sam_index1 = 0; /* sam1 index into probs matrix. */
  int pr_mat_sam_index2 = 0; /* sam2 index into probs matrix. */
  double numer_sum = 0.0;  /* Sums for each sample. */
  double* mags = 0; /* The magnitude of each 8 founder vector
                     * for each sample at each SNP. num_samples * num_snps. */
  double sum = 0.0;
  
  /* Zero out the K matrix. */
  for(sam1 = 0; sam1 < num_samples * num_samples; sam1++) {
	K[sam1] = 0.0;
  } //* for(sam1) */
  
  /* Calculate the magnitude of each vector. */
  mags = calloc(num_samples * num_snps, sizeof(double));
  for(snp = 0; snp < num_snps; snp++) {
    if(snp % 1000 == 0)	R_CheckUserInterrupt();
    pr_mat_snp_index = snp * pr_mat_slice;
    for(sam1 = 0; sam1 < num_samples; sam1++) {
      sum = 0.0;
      for(f1 = 0; f1 < num_founders; f1++) {
        pr_mat_sam_index1 = pr_mat_snp_index + f1 * num_samples + sam1;
	     sum += probs[pr_mat_sam_index1] * probs[pr_mat_sam_index1];
      } /* for(f1) */
	  mags[snp * num_samples + sam1] = sqrt(sum);
	} /* for(sam) */
  } /* for(snp)*/
  
  /* Calculate cosine of the angle between each sample at each SNP and sum them.
   * We go through the data in this order because the probs are ordered in 
   * memory as samples * founders * SNPs.*/
  for(snp = 0; snp < num_snps; snp++) { 
    if(snp % 1000 == 0)	{
	  R_CheckUserInterrupt();
	  Rprintf("SNP %d\n", snp);
    } /* if(snp % 1000 == 0) */
    pr_mat_snp_index = snp * pr_mat_slice;
    for(sam1 = 0; sam1 < num_samples; sam1++) {
      for(sam2 = sam1; sam2 < num_samples; sam2++) {
        /* Calculate the dot product of the two vectors. */
        numer_sum = 0.0;
        for(f1 = 0; f1 < num_founders; f1++) {
          pr_mat_sam_index1 = pr_mat_snp_index + f1 * num_samples + sam1;
          pr_mat_sam_index2 = pr_mat_snp_index + f1 * num_samples + sam2;
          numer_sum += probs[pr_mat_sam_index1] * probs[pr_mat_sam_index2];
        } /* for(f1) */
		  /* Divide the dot product by the product of the magnitudes. */
        K[sam1 * num_samples + sam2] += numer_sum / mags[snp * num_samples + sam1] / 
		                                mags[snp * num_samples + sam2];
      } /* for(sam2) */
    }  /* for(sam1) */
  } /* for(snp) */

  /* Take the mean of each value. */
  for(sam1 = 0; sam1 < num_samples; sam1++) {
    for(sam2 = sam1; sam2 < num_samples; sam2++) {
      K[sam1 * num_samples + sam2] = K[sam1 * num_samples + sam2] / (double)num_snps;
      K[sam2 * num_samples + sam1] = K[sam1 * num_samples + sam2];
    } /* for(sam2) */
  }  /* for(sam1) */

  free(mags);
  
} /* kinship() */


