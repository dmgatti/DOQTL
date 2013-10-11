/**********************************************************************
 * ri4hap.c
 *
 * functions to calculate the haplotype probabilities at generation
 * G1Fk in 4-way RIL by sibling mating
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
 * probability of AA haplotype on autosome at generation G1Fk in 
 * 4-way RIL by sibling mating
 * 
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate 
 * generations in the construction of recombinant inbred lines 
 **********************************************************************/
double ri4way_auto_hapAA(double r, int k)
{
  double result, s, rsq;

  rsq = r*r;
  s = sqrt(4.0*rsq-12.0*r+5.0);

  if(r==0.5) {
    if(k==1) result = 1.0/8.0;
    else result = 1.0/16.0;
  }
  else {
    result = (0.25)*(1.0/(1.0+6.0*r) + 
		     (6.0*rsq-7.0*r+3.0*r*s)/((1.0+6.0*r)*s)*
		     pow((1.0 - 2.0*r - s)/4.0, (double)k) - 
		     (6.0*rsq-7.0*r-3.0*r*s)/((1.0+6.0*r)*s)*
		     pow((1.0 - 2.0*r + s)/4.0, (double)k));
  }

  return(result);
}


/**********************************************************************
 * probability of AA haplotype on X chr in female at generation G1Fk in 
 * 4-way RIL by sibling mating
 * 
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate 
 * generations in the construction of recombinant inbred lines 
 **********************************************************************/
double ri4way_femX_hapAA(double r, int k)
{
  double result, t, rsq, rcube;

  rsq = r*r;
  rcube = r*rsq;
  t = sqrt(rsq-10.0*r+5.0);

  result = (0.5)*(2.0/(12.0*r+3.0) + 1.0/(3.0*r+3.0)*pow(-0.5,(double)k) - 
		  (4.0*rcube - t*(4.0*rsq+3.0*r)+3.0*rsq-5.0*r)/((8.0*rsq+10.0*r+2.0)*t)*
		  pow((1.0 - r + t)/4.0, (double)k) + 
		  (4.0*rcube + t*(4.0*rsq+3.0*r)+3.0*rsq-5.0*r)/((8.0*rsq+10.0*r+2.0)*t)*
		  pow((1.0 - r - t)/4.0,(double)k));

  return(result);
}


/**********************************************************************
 * probability of CC haplotype on X chr in female at generation G1Fk in 
 * 4-way RIL by sibling mating
 * 
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate 
 * generations in the construction of recombinant inbred lines 
 **********************************************************************/
double ri4way_femX_hapCC(double r, int k)
{
  double result, t, rsq;

  rsq = r*r;
  t = sqrt(rsq-10.0*r+5.0);

  result = 1.0/(12.0*r+3.0) - 1.0/(3.0*r+3.0)*pow(-0.5,(double)k) + 
    (9.0*rsq +5.0*r + r*t)/((8.0*rsq+10.0*r+2.0)*t)*pow((1.0 - r + t)/4.0,(double)k) - 
       (9.0*rsq +5.0*r - r*t)/((8.0*rsq+10.0*r+2.0)*t)*pow((1.0 - r - t)/4.0,(double)k);

  return(result);
}


/**********************************************************************
 * probability of AA haplotype on X chr in male at generation G1Fk in 
 * 4-way RIL by sibling mating
 * 
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate 
 * generations in the construction of recombinant inbred lines 
 **********************************************************************/
double ri4way_malX_hapAA(double r, int k)
{
  double result, t, rsq, rcube, rfourth;

  rsq = r*r;
  rcube = r*rsq;
  rfourth = rsq*rsq;
  t = sqrt(rsq-10.0*r+5.0);

  result = 1.0/(12.0*r+3.0) - 1.0/(3.0*r+3.0)*pow(-0.5, (double)k) + 
    (rcube - t*(8.0*rcube+rsq-3.0*r)-10.0*rsq+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5.0)/2 *
    pow((1.0 - r + t)/4.0,(double)k) + 
    +(rcube + t*(8.0*rcube+rsq-3.0*r)-10.0*rsq+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5.0)/2 *
    pow((1.0 - r - t)/4.0,(double)k);

  return(result);
}


/**********************************************************************
 * probability of CC haplotype on X chr in male at generation G1Fk in 
 * 4-way RIL by sibling mating
 * 
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate 
 * generations in the construction of recombinant inbred lines 
 **********************************************************************/
double ri4way_malX_hapCC(double r, int k)
{
  double result, t, rsq, rcube, rfourth;

  rsq = r*r;
  rcube = r*rsq;
  rfourth = rsq*rsq;
  t = sqrt(rsq-10.0*r+5.0);

  result = 1.0/(12.0*r+3.0) + 2.0/(3.0*r+3.0)*pow(-0.5, (double)k) + 
    (2.0*rfourth + t*(2.0*rcube-rsq+r)-19.0*rcube+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5) *
    pow((1.0 - r + t)/4,(double)k) + 
    (2.0*rfourth - t*(2.0*rcube-rsq+r)-19.0*rcube+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5) *
    pow((1.0 - r - t)/4,(double)k);

  return(result);
}


/* end of ri4hap.c */
