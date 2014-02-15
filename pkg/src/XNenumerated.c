//
//  XNenumerated.c
//
//  Created by Bill Engels on 12/16/13.
//  Copyright (c) 2013 Bill Engels. All rights reserved.
//
/*
 *     DESCRIPTION
 * Perform exact test for multinomial with full enumeration
 * The probabilities for each category in the multinomial are assumed to be known.
 * All possible outcomes of a multinomial experiment are examined with the aim of computing the
 * total probability of any outcome at least as distant from the expectation as observed.
 *
 * Three measures of distance from expectation are used:
 *   Log Likelihood Ratio (LLR)
 *   The probability of the outcome itself
 *   The standard chi square value
 *
 * A recursion method is used to find all possible outcomes. The three statistics are computed
 * in a distributive fashion for optimum performance.
 *
 *     EXAMPLE
 * Mendel's data from his 'dihybrid cross' : He observed 315, 108, 101,32 of the four seed types
 * in a cross from which he expected to see in the ratio 9:3:3:1. We want the probability of seeing
 * results _at least_ as deviant from the expected ratio as those observed among the 28,956,759 possible
 * outcomes when 556 objects are classified into four types. The probability of each outcome is given by 
 * the standard multinomial with the expected proportions used to compute the probabilities.
 *
*/


#ifdef NOT_READY_FOR_R    // This stuff can be removed in the R-only version
    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include "XNenumerated.h"

    #define R_pow_di pow  // The following defines convert R calls to normal C
    #define Calloc(x,y) calloc(x,sizeof(y))
    #define Free free
    #define Rprintf printf
    #define lgammafn lgamma
double factorial (unsigned k) {
    if(k > 170) return -1;
    if (k > 1) return (double) k * factorial(k-1);
    return 1;
}

double binomialCalc (unsigned n, unsigned k) {
	double coef = 1;
	double ndouble = n;
	double kdouble = k;
	unsigned i;
	if (k > n-k) {
		for (i = 1; i <= n-k; i++) coef *= (kdouble + i)/i;
	} else {
		for (i = 1; i <= k; i++) coef *= (ndouble-kdouble+i)/i;
	}
	return coef;
}

double multinomialProbability(int * k, double * p, int m) {
    double x = 1, pn;
    int i, j, nkount = 1;
    for (i = 0; i < m; i++) {
        pn = p[i];
        for (j = 0; j < k[i]; j++) {
            x *= pn * nkount / (j+1.);
            nkount++;
        }
    }
    return x;
}

//If R sends a negative number (-k) as.integer and C received it as unsigned, it is 4294967296-k in C.
void binomialCoeff (unsigned * ns, unsigned * ks, double * bc) {
    *bc = binomialCalc(*ns, *ks);
}

#else
    #include <R.h>
    #include <Rmath.h>
#endif



// Global variables.  These must be global for recursion.
static double *probs;
static double *expected;
static double *mKount;
static double pvalLLR, pvalProb, pvalChi;
static double multiobsLLR, multiobsProb, multiobsChi;
static double probSum, probRatio, probPerfect;
static double *flookup0, *flookup1;
static double *clookup0, *clookup1;
static int statType; // 1 for LLR, 2 for prob, 3 for chisquare

// histogram globals
static int histoBins = 0; // histoBins is only nonzero if a histogram is being made
static double * histo; // array with histogram data for LLR, prob, etc
static double statLeft, statSpan; //  scale statistics for histogram.

// All three statistics are computed, but we only want to use one of them for the histogram to save memory
static void handleHisto(double prob, double LLRstat, double chistat) {
    int hdex;
    double stat = 0;
    switch (statType) {
        case 1:
            stat = -2. * LLRstat; // convert LLR to D statistic
            break;
        case 2:
            stat = -2. * log(prob/probPerfect);
            break;
        case 3:
            stat = chistat;
            break;
        default:
            break;
    }
    hdex = (stat - statLeft)/statSpan;
    if ((hdex >= 0) && (hdex < histoBins)) {
        histo[hdex] += prob;
    }
}

// Use recursion to find all the possible multinomial outcomes while distributively computing their statistics and probabilities
static void multinom (int n, int k, double lprob, double LLR, double ChiStat)
{
    double LLRtemp, Chitemp, prob, exptemp = expected[k-1];
	if(k > 2)
	{
		multinom(n, k-1, lprob, LLR, ChiStat + exptemp);
		for (int i = 1; i <= n; i++ )
		{
			lprob = lprob  + log( (n-i+1) * probs[k-1]/(i * probs[0]));
            Chitemp = exptemp - i;
            Chitemp *= Chitemp;
            Chitemp /= exptemp;
			multinom(n-i, k-1, lprob, LLR +  i* log(probs[k-1]/i), ChiStat + Chitemp);
		}
	}
	else  // k = 2
	{
        prob = exp(lprob);
        if (n == 0) {
            if (LLR <= multiobsLLR) pvalLLR += prob;
            if (prob <= multiobsProb) pvalProb += prob;
            Chitemp = ChiStat + clookup0[0] + clookup1[0];
            if (Chitemp >= multiobsChi) pvalChi += prob;
            if(histoBins) handleHisto(prob, LLR, Chitemp);
            return;
        } else {
            for (long i = 0; i <= n; i++) {  // i is the second-last position, n-i is  in the last position
                if (i != 0) prob *= probRatio * (n + 1 - i)/i;
                LLRtemp = LLR +flookup1[i] + flookup0[n-i];
                Chitemp = ChiStat + clookup1[i] + clookup0[n-i];
                if (LLRtemp <= multiobsLLR) pvalLLR += prob;
                if (prob <= multiobsProb) pvalProb += prob;
                if (Chitemp >= multiobsChi) pvalChi += prob;
                if (histoBins) handleHisto(prob, LLRtemp, Chitemp);
            }
        }
	}
}

// called from R via .C
void exactMultinomialTest (int * obs,
                      double * expr,
                      int * nn,
                      int * statTypeR,
                      double * pLLR, // the LLR p-value
                      double * pProb, // the prob p-value
                      double * pChi, // the chi sq p-value
                      double * obsLLR, // the observed LLR
                      double * obsProb, // the observed prob
                      double * obsChiStat, // observed Chi Sq statistic
                      int * histoBinsR,
                      double * histoBounds,
                      double * histoData) {
    
    
    // get the total sample size, n
    unsigned n = 0;
    for (int i = 0; i < (*nn); i++) n += obs[i];
    
    // get memory for arrays
    probs = Calloc((*nn), double);
    expected = Calloc((*nn), double);
    unsigned maxlookup = n + 1;
    flookup0 = Calloc(maxlookup, double);
    flookup1 = Calloc(maxlookup, double);
    clookup0 = Calloc(maxlookup, double);
    clookup1 = Calloc(maxlookup, double);
    
    // Get probs and expecteds. The ones provided might not sum to one, so normalize them
    double exum = 0;
    for (int i = 0; i < (*nn); i++) exum += expr[i];
    for (int i = 0; i < (*nn); i++) {
        probs[i] = expr[i]/exum;
        expected[i] = probs[i] * n;
    }
    
    // Put observed values into global variables
    multiobsChi = *obsChiStat;
    multiobsProb = *obsProb;
    multiobsLLR = *obsLLR;
    
#ifdef NOT_READY_FOR_R
    // Get the observed probability (the 'dmultinom' function in R has no C interface
    multiobsProb = multinomialProbability(obs, probs, *nn);
    *obsProb = multiobsProb;

    // Get the observed LLR and chisq
    multiobsLLR = (double)n * log(n);
    for (int i = 0; i < *nn; i++) if(obs[i] > 0) multiobsLLR += obs[i] * log(probs[i]/obs[i]);
    *obsLLR = multiobsLLR;
    if (multiobsLLR > 0) {
        multiobsLLR = 1;   // If calculation is positive it means obs fits exp exactly, and roundoff error occurs
    }
#endif

    
    // Place the largest probability first to avoid underflow
    int bindex = 0;
    for (int i = 1; i < *nn; i++) {
        if (probs[bindex] < probs[i]) {
            bindex = i;
        }
        double probtemp = probs[0], expectedtemp = expected[0];
        probs[0] = probs[bindex];
        expected[0] = expected[bindex];
        probs[bindex] = probtemp;
        expected[bindex] = expectedtemp;
    }
    
    // find the best-fitting outcome's probability
    probPerfect = lgammafn(1. + n);
    int intexpi;
    for (int i = 0; i < *nn; i++) {
        intexpi = round(expected[i]);
        probPerfect += intexpi * log(probs[i]) - lgammafn(1. + intexpi);
    }
    probPerfect = exp(probPerfect);
    
    // Adjust the observed to avoid tests for floating equality
    double adj = 1.0000000001;
    multiobsProb *= adj;
    multiobsLLR /= adj;
    multiobsChi /= adj;
    
    
    // make lookup tables for faster execution
    flookup0[0] = 0;
    flookup1[0] = 0;
    
    for (int i = 1; i <= n; i++) {
        flookup0[i] = i * log(probs[0]/i);
        flookup1[i] = i * log(probs[1]/i);
    }
    for (int i = 0; i <= n; i++) {
        clookup0[i] = R_pow_di((expected[0] - i), 2)/expected[0];
        clookup1[i] = R_pow_di(expected[1] - i, 2)/expected[1];
    }
    
    if (*histoBinsR) // Prepare for histogram
    {
        histoBins = *histoBinsR;
        histo = histoData;
        for (int i = 0; i < histoBins; i++) histo[i] = 0.;
        statLeft = histoBounds[0];
        statSpan = (histoBounds[1] - statLeft)/(histoBins);
        if(statSpan == 0) histoBins = 0; // No histogram can be made
        statType = *statTypeR;
    }
    
    // Initialize some variables before call to multinom
    pvalLLR = 0; pvalProb = 0; pvalChi = 0;
    double kownt;
    mKount = &kownt;
    probSum = 0;
    probRatio = probs[1]/probs[0];
    int ncats = *nn;
    
    //********************************   Call the recursive routine to do all the work!
    
    multinom(n, ncats, log(probs[0]) * n, (double)n * log(n), 0);
    
    //********************************
    
    // Put the resulting three P-values into the return space
    *pLLR = pvalLLR;
    *pProb = pvalProb;
    *pChi = pvalChi;

    // Release memory
    Free(flookup0);
    Free(flookup1);
    Free(clookup0);
    Free(clookup1);
    Free(probs);
    Free(expected);
}
