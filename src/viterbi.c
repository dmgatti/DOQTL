/*******************************************************************************
 *  viterbi.c
 *  
 *  Created by Daniel Gatti on 4/11/11.
 *  Copyright 2011 The Jackson Laboratory. All rights reserved.
 *  int* dims: vector with three values: number of states, number of samples,
 *             number of SNPs.
 *  double* a: 3D matrix of transition probilities. num states * num states *
 *             num SNPs.
 *  double* b: 3D matrix of emission log-probabilities. num states * num samples
 *              * num SNPs.
 *  int* v_path: matrix with integer values corresponding to the states of
 *               the Viterbi path for each sample (num_samples * num_snps).
 *  double* v_prob: matrix with double values corresponding to the probability of
 *                  the Viterbi path for each sample (num_samples * num_snps).
 * R passes a pointer down for arrays.  We have to index into the 3D arrays by
 * hand.  SNPs are in the first slice.  Each slice has states in rows and 
 * samples in columns.  Since R uses column-ordering for matrices, indexing
 * through snps, then samples, then states should minimize cache-misses.
 * NOTE: The states are indexes from 0 to num_states-1.  In R, you will have to 
 *       add 1 to these values because R indexes from 1.
 */
#include <R.h>
#include <string.h>
#include "addlog.h"

void viterbi_from_r(int* dims, double* a, double* b, int* v_path, double* v_prob) {
  int snp = 0; /* index for SNPs */
  int sam = 0; /* index for samples */
  int st1 = 0; /* index for states */
  int st2 = 0; /* index for states */
  int num_states  = dims[0]; /* Total number of samples in prsmth. */
  int num_samples = dims[1];
  int num_snps    = dims[2];
  int outerindex = 0;  /* Index into 3D arrays. */
  int innerindex = 0;  /* Index into 3D arrays. */
  int slice_dim = num_states * num_samples;
  double max_prob = -DBL_MAX;  /* Maximum state for each sample. */
  int max_state   = 0;         /* Maximum state for each sample. */
  double init[num_states]; /* Initial probability. */
  double homo_init = log(1.0 / 64.0);
  double het_init  = log(1.0 / 32.0);
  double curr_pr = 0.0;        /* The current probability. */
  /* Probability of being in each state.  num_samples rows x num_snps columns. */
  double* pr = (double*)malloc(num_snps * num_samples * sizeof(double));
  /* Index to the state the produced the maximum probability at the 
   * previous state. num_samples rows x num_snps columns. */
  int* ptr = (int*)malloc(num_snps * num_samples * sizeof(int));
  
  memset(pr, 0, sizeof(pr));
  memset(ptr, 0, sizeof(ptr));
  
  Rprintf("A\n");
  
  /* Initialize probabilities. */
  for(st1 = num_states - 1; st1 >= 0; st1--) {
    init[st1] = het_init;
  } /* for(st1) */
  init[0] = homo_init;
  init[8] = homo_init;
  init[15] = homo_init;
  init[21] = homo_init;
  init[26] = homo_init;
  init[30] = homo_init;
  init[33] = homo_init;
  init[35] = homo_init;

  /* Initialize the starting values. */
  for(sam = 0; sam < num_samples; sam++) {
    /* outerindex moves to each sample at the first SNP. */
	outerindex = sam * num_states;
    for(st1 = 0; st1 < num_states; st1++) {
      /* innerindex moves to each state at a sample. */
	  innerindex = st1 + outerindex;
      pr[innerindex] = init[st1] + b[innerindex];
	} /* for(st1) */
  } /* for(sam) */
  
  Rprintf("B\n");

  /* Loop through each SNP, finding the most probable path to each state and
   * keeping an index to the state that we came from. */
  for(snp = 1; snp < num_snps; snp++) {
    for(sam = 0; sam < num_samples; sam++) {
	  /* outerindex moves us to the current sample at the current SNP. */
	  outerindex = sam * num_states + snp * slice_dim;
	  /* st1 is the state that we're going to. */
	  for(st1 = 0; st1 < num_states; st1++) {
	    /* innerindex moves us to the current state at the current sample. */
        innerindex = st1 + outerindex;
        max_prob = -DBL_MAX;
		max_state = 0;
		/* st2 is the state that we're coming from. */
		for(st2 = 0; st2 < num_states; st2++) {
		  curr_pr = pr[st2 + outerindex - slice_dim] + 
		            a[st2 + st1 * num_states] + b[innerindex];
		  if(curr_pr > max_prob) {
            max_prob  = curr_pr;
			max_state = st2;  /* Most probable source state. */
		  } /* if(curr_pr > max_prob) */
        } /* for(st2) */
		pr[innerindex]  = max_prob;
		ptr[innerindex] = max_state;
	  } /* for(st1) */
	} /* for(sam) */
  } /* for(snp) */

  Rprintf("C\n");

  /* Starting at the end, go back and find the most probable state path for
   * each sample.  Place the results in v_path. */
  snp = num_snps - 1;
  for(sam = 0; sam < num_samples; sam++) {
    /* outerindex moves us to the current sample at the current SNP. */
    outerindex = sam * num_states + snp * slice_dim;
    max_prob = -DBL_MAX;
    max_state = 0;
    /* Find the final state with the maximum probability. */
    for(st1 = 0; st1 < num_states; st1++) {
      if(pr[st1 + outerindex] > max_prob) {
        max_prob  = pr[st1 + outerindex];
	    max_state = st1;
      } /* if(pr[st1 + outerindex] > max_prob) */
    } /* for(st1) */
    /* Now outerindex takes us to the current sample and SNP. */
    outerindex = sam + snp * num_samples;
    v_prob[outerindex] = max_prob;
    v_path[outerindex] = max_state;
  } /* for(sam) */

  Rprintf("D\n");
  
  for(snp = num_snps - 2; snp >= 0; snp--) {
    for(sam = 0; sam < num_samples; sam++) {
      /* outerindex moves us to the current SNP and sample in v_prob/v_path. */
      outerindex = sam + snp * num_samples;
      /* innerindex moves us to the current SNP and sample in pr/ptr. */
      innerindex = sam + snp * num_samples + ptr[outerindex + num_samples];
      v_prob[outerindex] = pr[innerindex];
      v_path[outerindex] = ptr[innerindex];
    } /* for(sam) */
  } /* for(snp) */
  
  Rprintf("E\n");
  
  free(pr);
  free(ptr);

} /* viterbi_from_r() */

