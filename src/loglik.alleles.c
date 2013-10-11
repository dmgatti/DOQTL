/*******************************************************************************
 *  loglik.alleles.c
 *  
 *  Created by Daniel Gatti on Apr. 10, 2013.
 *  Copyright 2013 The Jackson Laboratory. All rights reserved.
 *  Calculate the log-likelihood of the model for the allele call case.
 *  int* dims: vector with four values: number of states, number of samples,
 *             number of SNPs, number of symbols.
 *  double* b: 3D array of emission probabilities. num_symbols x num_states x 
 *          num_snps.
 *  int* geno: 2D matrix of genotypes for each sample, coded as 0:3. 
 *             0 = AA, 1 = het, 2 = BB, 3 = nocall. num_samples x num_snps.
 *  double* prpred: 3D matrix of predictive log-probabilities. num states *
 *                  num samples * num SNPs.
 *  double* ll: Log-likelihood to return.
 *  
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "addlog.h"
void loglik_alleles(int* dims, double* b, int* geno, double* prpred, double* ll) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st  = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of states in prsmth. */
  int num_samples = dims[1]; /* Total number of samples in prsmth. */
  int num_snps    = dims[2]; /* Total number of SNPs in prsmth. */
  int num_symbols = dims[3]; /* Total number of symbols in prsmth. */
  int prpred_slice = num_states * num_samples; /* One SNP slice in prpred. */
  int b_slice      = num_states * num_symbols; /* One SNP slice in b. */
  int b_snp_index = 0;    /* Index into b array. */
  int geno_snp_index = 0; /* SNP index into geno array. */
  int geno_sam_index = 0; /* Sample index into geno array. */
  int pr_snp_index = 0;   /* SNP index into prpred array. */
  int pr_sam_index = 0;   /* Sample index into prpred array. */
  double sum = 0.0;       /* Accumulator for sum. */

  ll[0] = 0.0;
  for(snp = 0; snp < num_snps; snp++)
  {
    b_snp_index  = snp * b_slice;
	pr_snp_index = snp * prpred_slice;
	geno_snp_index = snp * num_samples;
    for(sam = 0; sam < num_samples; sam++) {
      sum = 0.0;
	  pr_sam_index = pr_snp_index + sam * num_states;
	  geno_sam_index = geno_snp_index + sam;
      for(st = 0; st < num_states; st++) {
        sum += exp(b[b_snp_index + st * num_symbols + geno[geno_sam_index]] + 
        		prpred[pr_sam_index + st]);
      } /* for(st) */
      ll[0] += log(sum);
    } /* for(sam) */
  } /* for(snp) */

} /*loglik_alleles() */
