/*******************************************************************************
 *  transition.c
 *  
 *  Functions to produce all transition matrix probabilities.  Recombination
 *  code written by Karl W. Broman, Aug. 23, 2011 (DOstep.c).
 *  
 *  Created by Daniel Gatti on 10/12/11.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with four values: number of states, number of samples,
 *             number of SNPs, number of samples to use.
 *  int* use: a vector of integer sample indices to use.  This allows some samples
 *            to be excluded from filtering (ie. founder samples).
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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "DOstep.h"

#include "DOstep.h"

/**********************************************************************
 * probability of recombinant haplotype on autosome at generation s 
 * of the diversity outcross population, where at generation 1 the 
 * mice are the result of random crosses among pre-CC mice with alpha_k[k]
 * denoting the proportion of those mice that are at generation k
 * NOTE: R passes everything down as a pointer, so the int values are 
 *       pointers and must be dereferenced.
 * Arguments:
 * double* r: vector of recombination fractions between SNPs.
 * int* s: DO outbreeding generation.
 * int* n_k: length of k and alpha_k
 * int* k: preCC inbreeding generation.
 * double* alpha_k: proportion of preCC progenitors at generation k.
 * int* num_recomb: length of recomb (and r) vector.
 * Returns:
 * double* recomb: with recombination probabilities for each interval in r.
 *                 NOTE: this must be pre-allocated in R.
 * From eqn 6 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
void DO_autosome_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k, 
		                     int* num_recomb, double* recomb)
{
  double hapAA = 0.0;
  int i = 0;
  int j = 0;

  for(i = 0; i < *num_recomb; i++) {
    /* Calculate probability of AA haplotype at generation s = 1 */
	hapAA = 0.0;
    for(j = 0; j < *n_k; j++) { 
      hapAA += (alpha_k[j] * ri4way_auto_hapAA(r[i], k[j] + 1) * (1.0 - r[i]) / 2.0);
    } /* for(j) */

    /* Later generations (using eqn (6)) */
    if(*s > 1)
      hapAA = 1.0 / 64.0 + pow(1.0 - r[i], (double)(*s - 1)) * (hapAA - 1.0 / 64.0);
  
    /* Probability of recombinant, assuming random order of initial crosses */
    recomb[i] = 1.0 - 8.0 * hapAA;
  } /* for(i) */

} /* DO_autosome_recomb_freq() */


/********************************************************************************
 * Probability of recombinant haplotype on female X chr at generation s 
 * of the diversity outcross population, where at generation s the 
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 * 
 * From eqn 9 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
void DO_femaleX_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k,
                            int* num_recomb, double* recomb)
{
  double f1, m1;
  double z, ws, ys;
  int i = 0;

  for(i = 0; i < *num_recomb; i ++ ) {
    if(*s == 1) 
      recomb[i] = DOrec_femX_s1(r[i], *n_k, k, alpha_k);
    else {
      z = sqrt((1.0 - r[i]) * (9.0 - r[i]));
      ws = pow((1.0 - r[i] + z) / 4.0, (double)(*s-1));
      ys = pow((1.0 - r[i] - z) / 4.0, (double)(*s-1));

      /* calculate probability of AA haplotype at generation s=1 */
      f1 = DOrec_femX_s1(r[i], *n_k, k, alpha_k);
      m1 = DOrec_malX_s1(r[i], *n_k, k, alpha_k);

      recomb[i] = (2.0 + (-64.0 * f1 * (1.0 - r[i]) - 128.0 * m1 + 3.0 - r[i]) / z *
    		      (ys - ws) - (1.0 - 64 * f1) * (ws + ys)) / 128.0;
    } /* else */
    recomb[i] = 1.0 - 8.0 * recomb[i];
  } /* for(i) */
  
} /* DO_femaleX_recomb_freq() */



/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s 
 * of the diversity outcross population, where at generation s the 
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 * 
 * From eqn 9 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
void DO_maleX_recomb_freq(double* r, int* s, int* n_k, int *k, double *alpha_k,
	                      int* num_recomb, double* recomb)
{
  double f1, m1;
  double z, ws, ys;
  int i = 0;

  for(i = 0; i < *num_recomb; i ++ ) {
    if(*s == 1) 
      recomb[i] = DOrec_malX_s1(r[i], *n_k, k, alpha_k);
    else {
      z = sqrt((1.0 - r[i]) * (9.0 - r[i]));
      ws = pow((1.0 - r[i] + z) / 4.0, (double)(*s - 1));
      ys = pow((1.0 - r[i] - z) / 4.0, (double)(*s - 1));

      /* calculate probability of AA haplotype at generation s=1 */
      f1 = DOrec_femX_s1(r[i], *n_k, k, alpha_k);
      m1 = DOrec_malX_s1(r[i], *n_k, k, alpha_k);

      recomb[i] = (2.0 + (64.0 * m1 - 256.0 * f1 + 3.0) * (1.0 - r[i]) / z * 
    		      (ys - ws) - (1.0 - 64 * m1) * (ws + ys)) / 128.0;
    } /* else */
 
    recomb[i] = 1.0 - 8.0 * recomb[i];
  } /* for(i) */
  
} /* DO_maleX_recomb_freq() */



