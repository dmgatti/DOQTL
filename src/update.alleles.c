/*******************************************************************************
 *  update_alleles.c
 *  
 *  Created by Daniel Gatti on Apr. 10, 2013.
 *  Copyright 2013 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with four values: number of states, number of samples,
 *             number of SNPs, number of symbols.
 *  int* geno: 2D matrix of genotypes for each sample, coded as 0:3. 
 *             0 = AA, 1 = het, 2 = BB, 3 = nocall. num_samples x num_snps.
 *  double* b: 3D array of log emission probabilities. num_symbols x 
 *             num_states x num_snps.
 *  double* pseudo: 3D array of log pseudocounts. num_symbols x num_states x 
 *          num_snps. 
 *  double* prsmth: 3D matrix of smoothed log-probabilities. num states *
 *                  num samples * num SNPs.
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through states, then samples, then SNPs should minimize cache-misses.
 * At each SNP, the variance is pooled among all states and samples.
 * We perform the calculations on a non-log scale.
 */
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "addlog.h"

void update_alleles(int* dims,int* geno, double* b, double* pseudo, 
                           double* prsmth) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st  = 0; /* index for states */
  int sym = 0; /* index for symbols */
  int num_states  = dims[0]; /* Total number of states in prsmth. */
  int num_samples = dims[1]; /* Total number of samples in prsmth. */
  int num_snps    = dims[2]; /* Total number of SNPs in prsmth. */
  int num_symbols = dims[3]; /* Total number of symbols in b. */
  int prsmth_slice = num_states * num_samples; /* One SNP slice in prsmth. */
  int b_slice = num_states * num_symbols;  /* One SNP slice in b. */
  int prsmth_snp_index = 0; /* SNP index into prsmth. */
  int prsmth_sam_index = 0; /* Sample index into prsmth. */
  int geno_snp_index = 0;   /* SNP index into geno. */
  int geno_value = 0;       /* Genotype value at current sample & SNP. */  
  int b_snp_index = 0;      /* SNP index into b array. */
  int b_index = 0;          /* Index into b array. */
  double b_sums[num_symbols]; /* Sum of symbols for each state. */
  double prsmth_sums = 0.0;   /* prsmth sums. */
  double pad = log(10.0);  /* The number of founder samples that we use to 
                              pad the emission probabilities. */

   /* Update the emission probabilities. 
    * We use all samples to update the emission probabilities. */
   for(snp = 0; snp < num_snps; snp++) {
      if(snp % 100 == 0) R_CheckUserInterrupt();

      /* Update the state means. */
      /* snp_index moves us to the current SNP in the mean arrays. */
 	  prsmth_snp_index = snp * prsmth_slice;
	  b_snp_index = snp * b_slice;
	  geno_snp_index = snp * num_samples;

	  /* Initialize the b_sums with the pseudocounts and zero out 
	   * the accumulators. */
      for(st = 0; st < num_states; st++) {
         b_index = b_snp_index + st * num_symbols;
         /* Pad the counts w/ the pseudocounts. */
         for(sym = 0; sym < num_symbols; sym++) {
            b_sums[sym] = pseudo[b_index + sym] + pad;
         } /* for(sym) */

         prsmth_sums = pad;

	     /* Use all samples to calculate the new emission probabilities. */
         for(sam = 0; sam < num_samples; sam++) {
        	prsmth_sam_index = st + sam * num_states + prsmth_snp_index;
			geno_value = geno[geno_snp_index + sam];
	        /* Update the sums for each symbol. */
            b_sums[geno_value] = addlog(b_sums[geno_value], 
			                     prsmth[prsmth_sam_index]);

            /* Add the current prsmth values to the sum for all samples at this
             * state. */
	        prsmth_sums = addlog(prsmth_sums, prsmth[prsmth_sam_index]);
         } /* for(sam) */
 
         /* Divide by the sum of the prsmth at this state. */
         for(sym = 0; sym < num_symbols; sym++) {
           b[b_index + sym] = b_sums[sym] - prsmth_sums;
         } /* for(sym) */

      } /* for(st) */
   } /* for(snp) */
} /* update_alleles() */

