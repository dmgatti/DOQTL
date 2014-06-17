/**********************************************************************
 * DOstep.c
 *
 * functions to calculate the step (aka transition) probabilities for
 * the diversity outcross population
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
 * transition probability for DO, autosome
 *
 * left = genotype at left locus
 * right = genotype at right locus
 * r = recombination fraction
 * s = generation of DO
 * 
 * alpha_k = proportion of preCC progenitors at generation k
 * n_k = length of k and alpha_k
 *
 * This calculates Pr(right | left)
 *
 * code left, right = 0..63 
 * with: first allele =  left / 8  (0..7)
 *       second allele = left % 8  (0..7)
 **********************************************************************/
double DOstep_auto(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k)
{
  int left1, left2, right1, right2;
  double recprob;

  /* pull out alleles for left and right loci */
  left1 = left / 8;
  left2 = left % 8;
  right1 = right / 8;
  right2 = right % 8;

  /* probability of recombinant haplotype */
  recprob = DOrec_auto(r, s, n_k, k, alpha_k);

  if(left1 == left2) { 
    if(right1 == right2) {
      if(left1 == right1) { /* AA -> AA */
        return( (1.0 - recprob)*(1.0 - recprob) );
      }
      else { /* AA -> BB */
        return( recprob*recprob/49.0 );
      }
    }
    else {
      if(left1 == right1 || left1 == right2) { /* AA -> AB */
        return( 2.0*recprob*(1.0-recprob)/7.0 );
      }
      else { /* AA -> BC */
        return( recprob * recprob * 2.0 / 49.0 );
      }
    }
  }
  else { /* AB */
    if(right1 == right2) {
      if(left1 == right1 || left2 == right1) { /* AB -> AA */ 
        return( recprob * (1.0 - recprob) / 7.0 );
      }
      else { /* AB -> CC */
        return( recprob * recprob / 49.0 );
      }
    }
    else {
      if((left1==right1 && left2==right2) ||
	 (left1==right2 && left2==right1)) { /* AB -> AB */
        return( recprob*recprob/49.0 + (1-recprob)*(1-recprob) );
      }
      else if(left1==right1 || left1==right2 ||
	      left2==right1 || left2==right2) { /* AB -> AC */
        return( recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0 );
      }
      else { /* AB -> CD */
        return( recprob*recprob*2.0/49.0 );
      }
    }
  }

} /* DOstep_auto() */

/**********************************************************************
 * transition probability for DO, female X chr
 *
 * code left, right = 0..63 
 * with: first allele =  left / 8 (0..7)
 *       second allele = left % 8 (0..7)
 **********************************************************************/
double DOstep_femX(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k)
{
  int left1, left2, right1, right2;
  double recprob;

  /* pull out alleles for left and right loci */
  left1 = left / 8;
  left2 = left % 8;
  right1 = right / 8;
  right2 = right % 8;

  /* probability of recombinant haplotype */
  recprob = DOrec_femX(r, s, n_k, k, alpha_k);

  if(left1 == left2) { 
    if(right1 == right2) {
      if(left1 == right1) { /* AA -> AA */
        return( (1.0 - recprob)*(1.0 - recprob) );
      }
      else { /* AA -> BB */
        return( recprob*recprob/49.0 );
      }
    }
    else {
      if(left1 == right1 || left1 == right2) { /* AA -> AB */
        return( 2.0*recprob*(1.0-recprob)/7.0 );
      }
      else { /* AA -> BC */
        return( recprob * recprob * 2.0 / 49.0 );
      }
    }
  }
  else { /* AB */
    if(right1 == right2) {
      if(left1 == right1 || left2 == right1) { /* AB -> AA */ 
        return( recprob * (1.0 - recprob) / 7.0 );
      }
      else { /* AB -> CC */
        return( recprob * recprob / 49.0 );
      }
    }
    else {
      if((left1==right1 && left2==right2) ||
	 (left1==right2 && left2==right1)) { /* AB -> AB */
        return( recprob*recprob/49.0 + (1-recprob)*(1-recprob) );
      }
      else if(left1==right1 || left1==right2 ||
	      left2==right1 || left2==right2) { /* AB -> AC */
        return( recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0 );
      }
      else { /* AB -> CD */
        return( recprob*recprob*2.0/49.0 );
      }
    }
  }
} /* DOstep_femX() */

/**********************************************************************
 * transition probability for DO, male X chr
 *
 * left, right coded 1, 2, ..., 8 [just one X chr in males]
 **********************************************************************/
double DOstep_malX(int left, int right, double r, int s, 
		   int n_k, int *k, double *alpha_k)
{
  double recprob;

  /* probability of recombinant haplotype */
  recprob = DOrec_malX(r, s, n_k, k, alpha_k); 

  if(left == right) return(1.0 - recprob);
  return(recprob/7.0);
} /* DOstep_malX() */

/* end of DOstep.c */
