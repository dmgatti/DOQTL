/**********************************************************************
 * DOstep_Rwrappers.c
 * 
 * wrappers for calling various functions related to calculate 
 * transitions probabilities for DO from R
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "DOstep.h"
#include "DOstep_Rwrappers.h"

void R_ri4way_auto_hapAA(double *r, int *k, double *result)
{
  *result = ri4way_auto_hapAA(*r, *k);
}

void R_ri4way_femX_hapAA(double *r, int *k, double *result)
{
  *result = ri4way_femX_hapAA(*r, *k);
}

void R_ri4way_femX_hapCC(double *r, int *k, double *result)
{
  *result = ri4way_femX_hapCC(*r, *k);
}

void R_ri4way_malX_hapAA(double *r, int *k, double *result)
{
  *result = ri4way_malX_hapAA(*r, *k);
}

void R_ri4way_malX_hapCC(double *r, int *k, double *result)
{
  *result = ri4way_malX_hapCC(*r, *k);
}

void R_DOrec_auto(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOrec_auto(*r, *s, *n_k, k, alpha_k);
}

void R_DOrec_femX(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOrec_femX(*r, *s, *n_k, k, alpha_k);
}

void R_DOrec_malX(double *r, int *s, int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOrec_malX(*r, *s, *n_k, k, alpha_k);
}


void R_DOstep_auto(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOstep_auto(*left, *right, *r, *s, *n_k, k, alpha_k);
}

void R_DOstep_femX(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOstep_femX(*left, *right, *r, *s, *n_k, k, alpha_k);
}

void R_DOstep_malX(int *left, int *right, double *r, int *s, 
		   int *n_k, int *k, double *alpha_k, double *result)
{
  *result = DOstep_malX(*left, *right, *r, *s, *n_k, k, alpha_k);
}


/* end of DOstep_Rwrappers.c */
