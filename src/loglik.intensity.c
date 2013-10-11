/*******************************************************************************
 *  loglik.intensity.c
 *  
 *  Created by Daniel Gatti on Apr. 10, 2013.
 *  Copyright 2013 The Jackson Laboratory. All rights reserved.
 *  Calculate the log-likelihood of the model for the intensity call case.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* b: 3D array of emission probabilities. nnum states *
 *             num samples * num SNPs.
 *  double* prpred: 3D matrix of predictive log-probabilities. num states *
 *                  num samples * num SNPs.
 *  double* ll: Log-likelihood to return.
 *  
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "addlog.h"
void loglik_intensity(int* dims, double* b, double* prpred, double* ll) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st  = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of states in prsmth. */
  int num_samples = dims[1]; /* Total number of samples in prsmth. */
  int num_snps    = dims[2]; /* Total number of SNPs in prsmth. */
  int slice = num_states * num_samples; /* One SNP slice in b or prpred. */
  int snp_index = 0;   /* SNP index into b & prpred arrays. */
  int sam_index = 0;   /* Sample index into b & prpred arrays. */
  int idx = 0;         /* Final array index for b & prpred. */
  double sum = 0.0;    /* Accumulator for sum. */

  ll[0] = -DBL_MAX;
  for(snp = 0; snp < num_snps; snp++)
  {
    snp_index  = snp * slice;
    for(sam = 0; sam < num_samples; sam++) {
	  sam_index = snp_index + sam * num_states;
      for(st = 0; st < num_states; st++) {
	    idx = sam_index + st;
        ll[0] += addlog(sum, b[idx] + prpred[idx]);
      } /* for(st) */
    } /* for(sam) */
  } /* for(snp) */

} /*loglik_intensity() */
