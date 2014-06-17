/**********************************************************************
 * DOrec.c
 *
 * functions to calculate the probability of a recombinant haplotype
 * in the diversity outcross population
 *
 * Karl W Broman
 * first written 23 Aug 2011
 * last modified 23 Aug 2011
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "DOstep.h"

/**********************************************************************
 * probability of recombinant haplotype on autosome at generation s 
 * of the diversity outcross population, where at generation 1 the 
 * mice are the result of random crosses among pre-CC mice with alpha_k[k]
 * denoting the proportion of those mice that are at generation k
 * 
 * From eqn 6 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_auto(double r, int s, int n_k, int *k, double *alpha_k)
{
  double hapAA;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  hapAA = 0.0;
  for(i=0; i<n_k; i++) 
    hapAA += (alpha_k[i] * ri4way_auto_hapAA(r, k[i]+1) * (1.0-r)/2.0);

  /* later generations (using eqn (6)) */
  if(s > 1)
    hapAA = 1.0/64.0 + pow(1.0-r, (double)(s-1)) * (hapAA - 1.0/64.0);
  
  /* probability of recombinant, assuming random order of initial crosses */
  return( 1.0 - 8.0*hapAA );
}

/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s 
 * of the diversity outcross population, where at generation 1 the 
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 * 
 * From eqn 7 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_femX_s1(double r, int n_k, int *k, double *alpha_k)
{
  double result;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  result = 0.0;
  for(i=0; i<n_k; i++) 
    result += (alpha_k[i] * (ri4way_femX_hapAA(r, k[i]+1) * (2.0-r) +
			     ri4way_femX_hapCC(r, k[i]+1) * (1.0-r)));

  return(result / 8.0);
}

/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s 
 * of the diversity outcross population, where at generation 1 the 
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 * 
 * From eqn 7 in Broman (2011) Haplotype probabilities in advanced 
 * intercross populations.  Technical report #223, Department of 
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_malX_s1(double r, int n_k, int *k, double *alpha_k)
{
  double result;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  result = 0.0;
  for(i=0; i<n_k; i++) 
    result += (alpha_k[i] * (ri4way_malX_hapAA(r, k[i]+1) * (2.0-r) +
			     ri4way_malX_hapCC(r, k[i]+1) * (1.0-r)));

  return(result / 8.0);
}



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
double DOrec_femX(double r, int s, int n_k, int *k, double *alpha_k)
{
  double result, f1, m1;
  double z, ws, ys;

  if(s==1) 
    result = DOrec_femX_s1(r, n_k, k, alpha_k);
  else {
    z = sqrt((1.0-r)*(9.0-r));
    ws = pow((1.0-r+z)/4.0, (double)(s-1));
    ys = pow((1.0-r-z)/4.0, (double)(s-1));

    /* calculate probability of AA haplotype at generation s=1 */
    f1 = DOrec_femX_s1(r, n_k, k, alpha_k);
    m1 = DOrec_malX_s1(r, n_k, k, alpha_k);

    result = (2.0 + (-64.0*f1*(1.0-r) - 128.0*m1 + 3.0 - r)/z * (ys - ws) -
	      (1.0 - 64*f1)*(ws+ys))/128.0;
  }

  return( 1.0 - 8.0*result );
}

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
double DOrec_malX(double r, int s, int n_k, int *k, double *alpha_k)
{
  double result, f1, m1;
  double z, ws, ys;

  if(s==1) 
    result = DOrec_malX_s1(r, n_k, k, alpha_k);
  else {
    z = sqrt((1.0-r)*(9.0-r));
    ws = pow((1.0-r+z)/4.0, (double)(s-1));
    ys = pow((1.0-r-z)/4.0, (double)(s-1));

    /* calculate probability of AA haplotype at generation s=1 */
    f1 = DOrec_femX_s1(r, n_k, k, alpha_k);
    m1 = DOrec_malX_s1(r, n_k, k, alpha_k);

    result = (2.0 + (64.0*m1 - 256.0*f1 + 3.0)*(1.0 - r)/z * (ys - ws) -
	      (1.0 - 64*m1)*(ws+ys))/128.0;
  }

  return( 1.0 - 8.0*result );
}

/* end of DOrec.c */
